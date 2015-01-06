/**********************************************************************

  2D3V PIC-MCC CODE '2D Arc-PIC'

  Copyright 2010-2014 CERN and Helsinki Institute of Physics.
  This software is distributed under the terms of the
  GNU General Public License version 3 (GPL Version 3),
  copied verbatim in the file LICENCE.md. In applying this
  license, CERN does not waive the privileges and immunities granted to it
  by virtue of its status as an Intergovernmental Organization
  or submit itself to any jurisdiction.

  Project website: http://arcpic.web.cern.ch/
  Developers: Helga Timko, Kyrre Sjobak

  arcbounds_original.cpp:
  Particle boundary conditions, same as in Helga's thesis.

***********************************************************************/

#include "arcbound_original.h"
#include "mydef.h"

#include "random.h"

#define   XTRN extern
#include "var.h"
#include "arrays1.h"
#undef XTRN

#include <cmath>

using namespace std;

ArcOriginal::ArcOriginal(vector<char*>& options) : ofile_arcboundsOriginalDat(NULL) {
  if (options.size() != 13) {
    cout << "Error in ArcOriginal(): Expected 13 options, got "
	 << options.size() << endl;
    exit(1);
  }
  
  //Interpret lines
  sscanf(options[0],  "%*[^:]%*[:] %lg", &(this->beta_tip));
  sscanf(options[1],  "%*[^:]%*[:] %lg", &(this->beta_f));
  sscanf(options[2],  "%*[^:]%*[:] %lg", &(this->j_melt));
  sscanf(options[3],  "%*[^:]%*[:] %lg", &(this->Remission));
  sscanf(options[4],  "%*[^:]%*[:] %lg", &(this->Remission_theor));
  sscanf(options[5],  "%*[^:]%*[:] %lg", &(this->SEY));
  sscanf(options[6],  "%*[^:]%*[:] %lg", &(this->r_Cu_e));

  char foo;
  sscanf(options[7],  "%*[^:]%*[:] %c",  &(foo));
  if      (foo == 'y') this->tipPi = true;
  else if (foo == 'n') this->tipPi = false;
  else {
    printf("Error in ArcOriginal::ArcOriginal(): tipPi has to be either 'y' or 'n' \n");
    exit(1);
  }
  
  sscanf(options[8],  "%*[^:]%*[:] %c",  &(foo));
  if      (foo == 'y') this->fracInjectStep = true;
  else if (foo == 'n') this->fracInjectStep = false;
  else {
    printf("Error in ArcOriginal::ArcOriginal(): fracInjectStep has to be either 'y' or 'n' \n");
    exit(1);
  }

  sscanf(options[9],  "%*[^:]%*[:] %lg", &(this->alpha_flat));
  
  sscanf(options[10],  "%*[^:]%*[:] %lg", &(this->heatspike_threshold));
  sscanf(options[11],  "%*[^:]%*[:] %lg", &(this->heatspike_yield));

  sscanf(options[12], "%*[^:]%*[:] %u",  &(this->file_timestep));

  //Initialization
  this->has_melted  = false;
  this->erosion     = 0;
  this->emitted     = 0;
  this->emitted_tip = 0;
}

ArcOriginal::~ArcOriginal() {
  if (ofile_arcboundsOriginalDat != NULL) fclose(ofile_arcboundsOriginalDat);
}

void ArcOriginal::print_par() const {
  printf( " - beta_tip:            %g \n",          this->beta_tip         );
  printf( " - beta_f:              %g \n",          this->beta_f           );
  printf( " - j_melt:              %g [A/cm^2]\n",  this->j_melt           );
  printf( " - Remission:           %g [dz]\n",      this->Remission        );
  printf( " - Remission_theor:     %g [dz]\n",      this->Remission_theor  );
  printf( " - SEY:                 %g \n",          this->SEY              );
  printf( " - r_Cu_e:              %g \n",          this->r_Cu_e           );
  printf( " - tipPi                %c \n",          this->tipPi ? 'y' : 'n');
  printf( " - fracInjectStep       %c \n", this->fracInjectStep ? 'y' : 'n');
  printf( " - alpha_flat:          %g \n",          this->alpha_flat       );
  printf( " - heatspike_threshold  %g [A/cm^2]\n",  this->heatspike_threshold);
  printf( "   = %g [particles/Ldb^2/ion_timestep]\n",
	  this->heatspike_threshold*(Ndb*dt_ion*Omega_pe) /
	  (6.7193e-12*n_ref*sqrt(T_ref)) );
  printf( "   = %g [particles/inner cell/ion_timestep]\n",
	  this->heatspike_threshold*(Ndb*dt_ion*Omega_pe) /
	  (6.7193e-12*n_ref*sqrt(T_ref)) * PI*SQU(dr) );
  printf( "   = %g [particles/2nd cell/ion_timestep]\n",
	  this->heatspike_threshold*(Ndb*dt_ion*Omega_pe) /
	  (6.7193e-12*n_ref*sqrt(T_ref)) * 3*PI*SQU(dr) );
  printf( " - heatspike_yield      %g [Cu/std. timestep (~1.77 fs)] = %i [superparticles / timestep]\n", 
	  this->heatspike_yield, (int) (this->heatspike_yield/N_sp * (Omega_pe/0.2 * sqrt(4e18/n_ref)) ) );
  printf( " - file_timestep        %u \n",          this->file_timestep    );
}

void ArcOriginal::init(unsigned int nr, double zmin, double zmax, double rmax) {
  ArcBounds::init(nr,zmin,zmax,rmax);
  v_inj_e = v_te*0.01;
  v_inj_i = vt_ions[2];

  this->sigma_heatspike = 1e6*nr;
}
void ArcOriginal::re_init(unsigned int nr, double zmin, double zmax, double rmax) {
  //ArcBounds::re_init(nr, zmin, zmax, rmax);
  ArcOriginal::init(nr, zmin, zmax, rmax); //ArcBounds::re_init just calls init()
}

void ArcOriginal::remove_i(Particle pa[], size_t &np, unsigned int sort) {
  if(sort != 1){
    if (np != 0) {
      printf("Error detected in ArcOriginal::remove_i(): %zu particles of sort=%u\n", np,sort);
      exit(1);
    }
    return;
  }

  size_t n_lost = 0;

  // check threshold of sputtering yield, only for ions, only at the cathode
  // 2D: flux will now depend on the area!! taken into account further below
  double threshold = heatspike_threshold*(Ndb*dt_ion*Omega_pe) / 
    (6.7193e-12*n_ref*sqrt(T_ref));  
  double current [nr];
  for (unsigned int i=0; i<nr; i++ ) current[i] = 0.;  

  //Store sputtering particles here, use different formulaes
  // to determine how much is actually sputtered depending on threshold
  static vector<Sput> sput_cathode_temp;

  //Total cathode yield (for sanity check)
  int Ycat_sum = 0;

  // SEY integer and fractional part
  int SEY_i = int ( SEY );
  double SEY_f = SEY - SEY_i;
  
  for (size_t n=0; n<np; n++ ) {
    //"Infinity"
    if ( pa[n].p.r >= rmax ) {
      removed_i[sort][2]++;
      removedIons.push_back(pa[n]);
      n_lost++;
      continue; 
    }
    //Cathode
    else if ( pa[n].p.z < zmin ) {
      removed_i[sort][0]++;
      current_i[sort][0] += 2*(int(pa[n].p.r))+1;
      current_cathode[ int(pa[n].p.r) ] += 1;
      removedIons.push_back(pa[n]);

      Sput newSput = calc_sput(pa[n], cs_ions[sort], current);
      if ( newSput.Y != 0 ) {
	sput_cathode_temp.push_back(newSput);
	Ycat_sum += newSput.Y;
      }

      // SEY = 0.5 = constant, only from ions hitting the cathode
      // SEY with registering r-coordinates
      // Set reasonable threshold for incident ion energy (e.g. 100 eV)
      if ( (SQU(pa[n].p.vz) + SQU(pa[n].p.vr) + SQU(pa[n].p.vt))*T_ref/(2*SQU(cs_ions[sort])) > 100 ) {
	if ( RAND <= SEY_f ) {
	  Sput foo;
	  foo.r = pa[n].p.r;
	  foo.Y = SEY_i+1;
	  sput_cathode_SEY.push_back(foo);
	}
	else if ( SEY_i > 0 ) {
	  Sput foo;
	  foo.r = pa[n].p.r;
	  foo.Y = SEY_i;
	  sput_cathode_SEY.push_back(foo);
	}
      }
      n_lost++;
      continue; 
    }
    //Anode
    else if ( pa[n].p.z >= zmax ) {	
      removed_i[sort][1]++;
      current_i[sort][1] += 2*(int(pa[n].p.r))+1;
      current_anode[ int(pa[n].p.r) ] -= 1;      
      removedIons.push_back(pa[n]);
      
      Sput newSput = calc_sput(pa[n], cs_ions[sort], NULL);
      if ( newSput.Y != 0 ) sput_anode.push_back(newSput);
      
      n_lost++;
      continue; 
    }
    //Implicit else: keep this particle
    pa[n-n_lost].p = pa[n].p;
  }
  np -= n_lost; 
  
  // Enhanced yield?
  bool check_enh = false;
  double r, sigma(1e6*nr); //Crash if sigma not initialized
  for (unsigned int i=0; i<nr; i++ ) {
    // Rescale current with area of cell / area Ldb^2
    current[i] /= PI * (2*i + 1) * SQU(dr);
    
    if ( current[i] >= threshold ) {
      check_enh = true;
      sigma = i+1;
    }
  }
  sigma_heatspike = sigma; //store for output
  
  if ( check_enh == false ) {
    // Yamamura-Tawara fitting for Cu -> Cu
    for (size_t j=0; j<sput_cathode_temp.size(); j++ ) {
      sput_cathode.push_back(sput_cathode_temp[j]);
    }
  }
  else if ( check_enh == 1 ) {
    // Enhanced yield from MD
    // Generate Gaussian distribution for r-coord's 
    int num_sputtered = heatspike_yield / N_sp * (Omega_pe/0.2 * sqrt(4e18/n_ref));
    for (int j=0; j<num_sputtered; j++ ) { // 1000/SN
      do { r = sigma * sqrt(-2.*log(RAND+1.e-20)); } while ( r >= nr );
      Sput foo;
      foo.r = r;
      foo.Y = 1;
      sput_cathode.push_back(foo);
    }
    
    // Print control message
    if ( Ycat_sum > 1000 ) { // 1000/SN
      printf("*** UNDERESTIMATED *** Yamamura sputtering yield Ycat_sum ( %d ) > enhanced sputtering yield ( %d )! \n", Ycat_sum, 1000 );
      fflush( stdout );
    } 
  }

  sput_cathode_temp.clear();
}
void ArcOriginal::remove_n(Particle pa[], size_t &np){
  size_t n_lost = 0; 
  
  for (size_t n=0; n<np; n++ ) {
    //"Infinity"
    if ( pa[n].p.r >= rmax ) {
      removed_n[2]++;
      removedNeutrals.push_back(pa[n]);
      n_lost++;
      continue; 
    }
    //Cathode
    else if ( pa[n].p.z < zmin ) {
      removed_n[0]++;
      current_n[0] += 2*(int(pa[n].p.r))+1;
      removedNeutrals.push_back(pa[n]);      

      Sput newSput = calc_sput(pa[n], cs_ions[NSpecies-1], NULL);
      if ( newSput.Y != 0 ) sput_cathode.push_back(newSput);

      n_lost++; 
      continue; 
    }
    //Anode
    else if ( pa[n].p.z >= zmax ) {
      removed_n[1]++;
      current_n[1] += 2*(int(pa[n].p.r))+1;
      removedNeutrals.push_back(pa[n]);

      Sput newSput = calc_sput(pa[n], cs_ions[NSpecies-1], NULL);
      if ( newSput.Y != 0 ) sput_anode.push_back(newSput);

      n_lost++; 
      continue; 
    }
    //Implicit else: keep this particle
    pa[n-n_lost].p = pa[n].p; 
  }
  np -= n_lost;  
}

Sput ArcOriginal::calc_sput(const Particle& pa, const double cs, double* current_enhancedY) {
  double nrg = SQU(pa.p.vz) + SQU(pa.p.vr) + SQU(pa.p.vt);
  nrg *= T_ref/(2*SQU(cs));
  if ( nrg > 23.383 ) { //in eV
    // Register fluxes going through each cell for enhanced Y
    if ( current_enhancedY != NULL ) {
      current_enhancedY[size_t(pa.p.r)] += 1.; // to be rescaled!
    }
    
    // Reduced energy
    double eps = (4.45394e-06)*nrg;
    
    // Nuclear stopping cross section
    double sn = 8205*3.441*sqrt(eps)*log(eps + 2.718) /
      (1 + 6.355*sqrt(eps) + eps*(6.882*sqrt(eps) - 1.708));
    
    // Sputtering yield 
    double p = 0.042*(0.2525/3.49) * (sn/(1 + 0.010124*0.1573*pow(eps,0.3))) *
      pow(1 - sqrt(23.383/nrg),2.5);
    
    // Fractional part of p handled probabilistically
    double Y = int ( p );
    p -= Y;
    if ( RAND <= p ) Y++;
    
    // Register r-coordinates of bombarding particles
    if ( Y > 0 ) {
      Sput ret;
      ret.r = pa.p.r;
      ret.Y = Y;
      return ret;
    }
  }
  
  //No sputtering
  Sput ret;
  ret.r = 0.0;
  ret.Y = 0;
  return ret;
}

void ArcOriginal::inject_e(Particle pa[], size_t &np, double const Ez[]) {
  //Field (dimless & GV/m inl. beta) and emitted currents(tip & flat)
  double field, Eloc, j;
  double jFN[nr];
  
  // FN AT FIELD EMITTER
  // Fowler-Nordheim with Wang-Loew approximation
  // W-L: v(y) ~ 0.956 - 1.062*(3.7947e-5)^2*Eloc/(4.5*4.5)
  // beta=dynamic, work fct=4.5eV, j in A/cm^2, E in V/m
  
  field = Ez[0];
  // rescale the field to GV/m
  Eloc = - 2.69036254e-10*dz/SQU(Omega_pe)*sqrt(T_ref*n_ref)*field*beta_tip; 
  this->E_tip_loc = Eloc; //For output

  if (field < 0.) {  
    // Protect against numerical fluctuations
#warning unphysical capping of emitted current
    if ( Eloc > 12. ) Eloc = 12.;
    
    j = 4.7133e9 * SQU(Eloc) * exp(-62.338/Eloc); // in A/cm^2
    j *= (tipPi ? PI : 1.0)*(Ndb*e2inj_step*Omega_pe*SQU(Remission_theor*dz))/(6.7193e-12*n_ref*sqrt(T_ref)); //dimless, 2D: rem^2
  }
  else { 
    j=0.; 
  }
  
  // FN OUTSIDE THE FIELD EMITTER
  for (unsigned int jj=0; jj<nr; jj++ ) {
    field = Ez[(jj+1)*NZ];
    if (field < 0.) {
      // B=B_f, rescale the field to GV/m
      Eloc = - 2.69036254e-10*dz/SQU(Omega_pe)*sqrt(T_ref*n_ref)*field*beta_f; 
      if ( Eloc > 12. ) Eloc = 12.;  
      jFN[jj] = 4.7133e9 * SQU(Eloc) * exp(-62.338/Eloc); // in A/cm^2

      jFN[jj] *= this->alpha_flat;

      //dimless currents
      if ( Remission <= jj )
	jFN[jj] *= (Ndb*e2inj_step*Omega_pe*PI*SQU(dz)*(2*jj+1))/(6.7193e-12*n_ref*sqrt(T_ref));
      else if ( (Remission > jj) && (Remission < jj+1) )
	jFN[jj] *= Ndb*e2inj_step*Omega_pe*PI*SQU(dz)*( SQU(jj+1)-SQU(Remission) )/(6.7193e-12*n_ref*sqrt(T_ref));
      else 
	jFN[jj] = 0.;
    }
    else { 
      jFN[jj]=0.; 
    }
  }
  
  //Done calculating how much should be injected: Now inject (only cathode) !

  double r1, r2; // temp variables from Gaussian RNG
  
  // SEY: inject with incident r-coordinates, empty SEY variables!
  size_t totY = 0;
  for (size_t k=0; k < sput_cathode_SEY.size(); k++ ) {
    for (int i=0; i < sput_cathode_SEY[k].Y; i++ ) {
      if ( ( np + totY ) >= NPART ) {
	printf("Error in ArcOriginal::inject_e(): Particle array overflow (SEY)\n");
	exit(1);
      }
      
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); }
      while( r1 > 5. );
      
      // Gaussian scheme
      r2 = RAND * TWOPI;
      pa[np+totY].p.vr = r1*cos(r2)*v_inj_e;
      pa[np+totY].p.vt = r1*sin(r2)*v_inj_e;
      
      // Gaussian scheme
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } 
      while( r1 > 5. );
      pa[np+totY].p.vz = r1*v_inj_e; 
      
      pa[np+totY].p.z = zmin + pa[np+totY].p.vz;
      pa[np+totY].p.r = sput_cathode_SEY[k].r + pa[np+totY].p.vr;
      if ( pa[np+totY].p.r < 0 ) pa[np+totY].p.r = 1.e-20;
      else if (pa[np+totY].p.r > nr) pa[np+totY].p.r = 2*nr - pa[np+totY].p.r;

      pa[np+totY].p.m = 1;
      
      current_e[0] += 2*( int(pa[np+totY].p.r) ) + 1;
      current_cathode[ int(pa[np+totY].p.r) ] += 1;

      totY++;
    }
  }
  np += totY;
  injected_e[0] += totY;
  sput_cathode_SEY.clear();

  // FN INJECTION FROM FIELD EMITTER
  // FE: inject with flat distribution over the emission radius
  //Inject 'tmp' num particles
  size_t tmp = size_t(j);
  j -= tmp;
  if ( RAND <= j ) tmp++;
  
  if ( ( np + tmp ) >= NPART) {
    printf("Error in ArcOriginal::inject_e(): Particle array overflow (FN at field emitter)\n");
    exit(1);
  }
  
  for  (size_t k=0; k<tmp; k++ ) {     
    if (fracInjectStep) {
      //Velocity
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
      r2 = RAND * TWOPI;
      pa[np+k].p.vr = r1*cos(r2)*v_inj_e;
      pa[np+k].p.vt = r1*sin(r2)*v_inj_e;
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
      pa[np+k].p.vz = r1*v_inj_e;
      
      //Position
      // uniform random position on disc
      pa[np+k].p.r = Remission*sqrt(RAND);
      pa[np+k].p.z = zmin;
      
      //Fractional timestep push (euler method)
      double Rp = RAND; //1-R; how (fractionally) far into the timestep
                        // are we at z=0?
      pa[np+k].p.r  += Rp*e2inj_step*pa[np+k].p.vr;
      pa[np+k].p.z  += Rp*e2inj_step*pa[np+k].p.vz;
      //No magnetic field, E = Ez on surface
      pa[np+k].p.vz -= (Rp*e2inj_step-0.5)*2*Ez[0];
      
      if ( pa[np+k].p.r < 0 ) pa[np+k].p.r = -pa[np+k].p.r; //Reflect on axis
      else if (pa[np+k].p.r > nr) pa[np+k].p.r = 2*nr - pa[np+k].p.r;
    }
    else {
      // Gaussian scheme
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
      r2 = RAND * TWOPI;
      pa[np+k].p.vr = r1*cos(r2)*v_inj_e;
      pa[np+k].p.vt = r1*sin(r2)*v_inj_e;
      
      // Gaussian scheme
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
      pa[np+k].p.vz = r1*v_inj_e; 
      
      pa[np+k].p.z = zmin + pa[np+k].p.vz;
      pa[np+k].p.r = Remission * sqrt(RAND) + pa[np+k].p.vr;
      
      if ( pa[np+k].p.r < 0 ) pa[np+k].p.r = 1.e-20;
      else if (pa[np+k].p.r > nr) pa[np+k].p.r = 2*nr - pa[np+k].p.r;
    }
    pa[np+k].p.m = 1;
    current_e[0] += 2*( int(pa[np+k].p.r) ) + 1;
    current_cathode[ int(pa[np+k].p.r) ] += 1;
  }
  np += tmp;

  injected_e[0] += tmp;  
  emitted       += tmp;
  emitted_tip   += tmp;
  
  // FN INJECTION OUTSIDE THE FIELD EMITTER
  for (unsigned int jj=0; jj<nr; jj++ ) {
    //Inject 'tmp' num. particles
    if ( jFN[jj] > 0.0 ) {
      tmp = (size_t) (jFN[jj]); 
      jFN[jj] -= tmp;
      if ( RAND <= jFN[jj] ) tmp++;
      
      if ( ( np + tmp ) >= NPART) {
	printf("Error in ArcOriginal::inject_e(): Particle array overflow (FN outside field emitter)\n");
	exit(1);
      }
      
      for  (size_t k=0; k<tmp; k++ ) { 
	if (fracInjectStep) {
	  //Velocity
	  do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
	  r2 = RAND * TWOPI;
	  pa[np+k].p.vr = r1*cos(r2)*v_inj_e;
	  pa[np+k].p.vt = r1*sin(r2)*v_inj_e;
	  do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
	  pa[np+k].p.vz = r1*v_inj_e;
	  
	  //Position
	  if ( Remission <= jj )
	    pa[np+k].p.r = jj + RAND;
	  else if ( (Remission > jj) && (Remission < jj+1) )
	    pa[np+k].p.r = Remission + (jj + 1 - Remission) * RAND;

	  pa[np+k].p.z = zmin;
	  
	  //Fractional timestep push (euler method)
	  double Rp = RAND; //1-R; how (fractionally) far into the timestep
	  // are we at z=0?
	  pa[np+k].p.r  += Rp*e2inj_step*pa[np+k].p.vr;
	  pa[np+k].p.z  += Rp*e2inj_step*pa[np+k].p.vz;
	  //No magnetic field, E = Ez on surface
	  pa[np+k].p.vz -= (Rp*e2inj_step-0.5)*2*Ez[(jj+1)*NZ];
	  
	  if ( pa[np+k].p.r < 0 ) pa[np+k].p.r = -pa[np+k].p.r; //Reflect on axis
	  else if (pa[np+k].p.r > nr) pa[np+k].p.r = 2*nr - pa[np+k].p.r;

	}
	else {	
	  // Gaussian scheme
	  do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
	  r2 = RAND * TWOPI;
	  pa[np+k].p.vr = r1*cos(r2)*v_inj_e; 
	  pa[np+k].p.vt = r1*sin(r2)*v_inj_e; 
	  
	  // Gaussian scheme	
	  do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
	  pa[np+k].p.vz = r1*v_inj_e;
	  
	  pa[np+k].p.z = zmin + pa[np+k].p.vz;
	  if ( Remission <= jj )
	    pa[np+k].p.r = jj + RAND + pa[np+k].p.vr;
	  else if ( (Remission > jj) && (Remission < jj+1) )
	    pa[np+k].p.r = Remission + (jj + 1 - Remission) * RAND + pa[np+k].p.vr;
	  
	  if ( pa[np+k].p.r < 0 ) pa[np+k].p.r = 1.e-20;
	  else if (pa[np+k].p.r > nr) pa[np+k].p.r = 2*nr - pa[np+k].p.r;
	}
	pa[np+k].p.m = 1;
	current_e[0] += 2*( int(pa[np+k].p.r) ) + 1;
	current_cathode[ int(pa[np+k].p.r) ] += 1;
      }
      np += tmp;

      injected_e[0] += tmp;
#warning Every FN electron counts for evaporation, but only injecting Cu close to tip
      emitted += tmp;
    }
  } 
}
void ArcOriginal::inject_n(Particle pa[], size_t &np, double const Ez[]) {
  // evaporation of neutrals
  double tmp = r_Cu_e * emitted;
  int n2inject_evap = int ( tmp );
  tmp -= n2inject_evap;
  if ( RAND <= tmp ) n2inject_evap++;
  emitted = 0;

  // Cathode, neutral evaporation
  double r1, r2;
  size_t m = n2inject_evap;
  if ( ( np + m ) >= NPART) {
    printf("Error in ArcOriginal::inject_n(): Particle array overflow (evaporation)\n");
    exit(1);
  }
  for ( size_t k=0; k<m; k++ ) {
    // Gaussian scheme
    do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
    r2 = RAND * TWOPI;
    pa[np+k].p.vr = r1*cos(r2)*v_inj_i; // do not suppress
    pa[np+k].p.vt = r1*sin(r2)*v_inj_i; // do not suppress
    
    // Gaussian scheme
    do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
    pa[np+k].p.vz = r1*v_inj_i;
    
    pa[np+k].p.z = zmin + pa[np+k].p.vz;
    pa[np+k].p.r = Remission * sqrt(RAND) + pa[np+k].p.vr;
    
    if ( pa[np+k].p.r < 0 ) pa[np+k].p.r = 1.e-20;
    else if (pa[np+k].p.r > nr) pa[np+k].p.r = 2*nr - pa[np+k].p.r;

    pa[np+k].p.m = 1; //SN;
    
  }
  np += m;
  erosion += n2inject_evap; //Erosion only due to neutral evap
  injected_n[0] += n2inject_evap;
  
  // Cathode, sputtering
  inject_sput(pa, np, sput_cathode, true);
  // Anode, sputtering only
  inject_sput(pa, np, sput_anode, false);
}
void ArcOriginal::inject_sput(Particle pa[], size_t &np, std::vector<Sput> &sput, bool isCathode) {
  
  double r1, r2;

  size_t m = sput.size();
  size_t tot = 0;
  for ( size_t k=0; k<m; k++ ) {
    for ( int i=0; i<sput[k].Y; i++ ) {
      if ( ( np + tot ) >= NPART) {
	printf("Error in ArcOriginal::inject_n(): Particle array overflow (sputtering)\n");
	exit(1);
      }
      
      // Gaussian scheme      
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
      r2 = RAND * TWOPI;
      pa[np+tot].p.vr = r1*cos(r2)*v_inj_i;
      pa[np+tot].p.vt = r1*sin(r2)*v_inj_i; 
      
      // Gaussian scheme
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
      pa[np+tot].p.vz = isCathode ? r1*v_inj_i : -r1*v_inj_i;
      
      pa[np+tot].p.z = (isCathode ? zmin : zmax) + pa[np+tot].p.vz;
      pa[np+tot].p.r = sput[k].r + pa[np+tot].p.vr; 
      
      if ( pa[np+tot].p.r < 0 ) pa[np+tot].p.r = 1.e-20;
      else if (pa[np+tot].p.r > nr) pa[np+tot].p.r = 2*nr - pa[np+tot].p.r;

      pa[np+tot].p.m = 1; //SN;
      
      // Sum up total yield for outputting
      tot++;
    }
  }
  np += tot;
  injected_n[isCathode ? 0 : 1] += tot;
  sput.clear();
}

void ArcOriginal::timestep(unsigned int nstep, bool isOutputTimestep) {
  //2D: use Rem (=emission radius) instead of 20nm
  if ( !has_melted ) {
    beta_tip -= erosion * n_ref / ( Ndb * pow(Remission_theor*dz,3)) * 3.824131e-24; 
    this->erosion = 0;
    if (beta_tip < beta_f) { beta_tip = beta_f; }
  }

  //Check for melting at same timesteps as Helga did (when doing output)
  // Melting of the emitter due to Ohmic heating by FE current, was .5e8
  // !!  This uses a ton of global variables (nsteps, nav_start, nav_time, dependency on dt_out,...)
  // !!  Should be *properly* decoupled from output!
  if ( nsteps == nav_start+nav_time &&
       !has_melted &&
       emitted_tip > (j_melt*Ndb*dt_out*SQU(Remission_theor*dz)*(tipPi ? PI : 1.0) / 
		      6.7193e-12/n_ref/sqrt(T_ref)) ) {
    beta_tip = beta_f;
    has_melted = true;   
  }
  
  //Write file and reset counting arrays
  ArcBounds::timestep(nstep, isOutputTimestep);
  emitted_tip = 0;
}

void ArcOriginal::writeFile_arcboundsOriginalDat(unsigned int nstep) {
  if (ofile_arcboundsOriginalDat == NULL) initFile_arcboundsOriginalDat();
  if ( nstep % file_timestep == 0) {
    fprintf(ofile_arcboundsOriginalDat, "%11d %8f               %c %27d %15f %19f\n", nstep, beta_tip, (has_melted? 'y' : 'n'), emitted_tip, E_tip_loc, sigma_heatspike);
    fflush(ofile_arcboundsOriginalDat);
  }
}
void ArcOriginal::initFile_arcboundsOriginalDat() {
  ofile_arcboundsOriginalDat = fopen("arcbounds_original.dat", "w");
  fprintf(ofile_arcboundsOriginalDat, "## timestep  beta_tip has_melted[y/n] emitted_tip[superparticles] E_tip_loc[GV/m] sigma_heatspike[dz]\n");
  fflush(ofile_arcboundsOriginalDat);
}
