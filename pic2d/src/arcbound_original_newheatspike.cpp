/**********************************************************************

  2D3V PIC-MCC CODE '2D Arc-PIC'

  Copyright 2010-2015 CERN and Helsinki Institute of Physics.
  This software is distributed under the terms of the
  GNU General Public License version 3 (GPL Version 3),
  copied verbatim in the file LICENCE.md. In applying this
  license, CERN does not waive the privileges and immunities granted to it
  by virtue of its status as an Intergovernmental Organization
  or submit itself to any jurisdiction.

  Project website: http://arcpic.web.cern.ch/
  Developers: Helga Timko, Kyrre Sjobak

  arcbounds_original_newheatspike.cpp:
  Particle boundary conditions, same as in Helga's thesis,
  except for the heat spike sputtering part which has been updated.

***********************************************************************/

#include "arcbound_original_newheatspike.h"
#include "mydef.h"

#include "random.h"

#define   XTRN extern
#include "var.h"
#include "arrays1.h"
#undef XTRN

#include <cmath>

using namespace std;

ArcOriginalNewHS::ArcOriginalNewHS(vector<char*>& options) : ofile_arcboundsOriginalDat(NULL) {

  if (options.size() != 17) {
    cout << "Error in ArcOriginalNewHS(): Expected 16 options, got "
	 << options.size() << endl;
    exit(1);
  }

  //Interpret lines
  char foo;

  sscanf(options[0],  "%*[^:]%*[:] %lg", &(this->beta_tip));
  sscanf(options[1],  "%*[^:]%*[:] %lg", &(this->beta_f));
  sscanf(options[2],  "%*[^:]%*[:] %lg", &(this->j_melt));

  sscanf(options[3],  "%*[^:]%*[:] %c",  &(foo));
  if      (foo == 'y') this->do_erosion = true;
  else if (foo == 'n') this->do_erosion = false;
  else {
    printf("Error in ArcOriginalNewHS::ArcOriginalNewHS(): do_erosion has to be either 'y' or 'n' \n");
    exit(1);
  }  
  if (this->do_erosion and (beta_tip < beta_f) ) {
    printf("Error in ArcOriginalNewHS::ArcOriginalNewHS(): when do_erosion = true, must have beta_tip > beta_flat\n"); 
    printf(" (else beta_tip gets set to beta_flat when doing the erosion)\n");
    exit(1);
  }

  sscanf(options[4],  "%*[^:]%*[:] %lg", &(this->Remission));
  sscanf(options[5],  "%*[^:]%*[:] %lg", &(this->Remission_theor));
  sscanf(options[6],  "%*[^:]%*[:] %lg", &(this->Rborder));

  sscanf(options[7],  "%*[^:]%*[:] %lg", &(this->SEY));
  sscanf(options[8],  "%*[^:]%*[:] %lg", &(this->r_Cu_e));
  sscanf(options[9],  "%*[^:]%*[:] %lg", &(this->r_Cu_e_flat));


  sscanf(options[10],  "%*[^:]%*[:] %c",  &(foo));
  if      (foo == 'y') this->evap_flat_center = true;
  else if (foo == 'n') this->evap_flat_center = false;
  else {
    printf("Error in ArcOriginalNewHS::ArcOriginalNewHS(): evap_flat_center has to be either 'y' or 'n' \n");
    exit(1);
  }
  
  sscanf(options[11],  "%*[^:]%*[:] %c",  &(foo));
  if      (foo == 'y') this->fracInjectStep = true;
  else if (foo == 'n') this->fracInjectStep = false;
  else {
    printf("Error in ArcOriginalNewHS::ArcOriginalNewHS(): fracInjectStep has to be either 'y' or 'n' \n");
    exit(1);
  }

  sscanf(options[12],  "%*[^:]%*[:] %lg", &(this->alpha_flat));

  sscanf(options[13],  "%*[^:]%*[:] %lg", &(this->heatspike_threshold));
  sscanf(options[14],  "%*[^:]%*[:] %lg", &(this->heatspike_yield_ratio));

  sscanf(options[15],  "%*[^:]%*[:] %c",  &(foo));
  if      (foo == 'a' or foo == 'b' or foo == 'c') {
    this->heatspike_model = foo;
  }
  else {
    printf("Error in ArcOriginalNewHS::ArcOriginalNewHS(): heatspike_model has to be either 'a', 'b' or 'c' \n");
    exit(1);
  }

  sscanf(options[16], "%*[^:]%*[:] %u",  &(this->file_timestep));

}

ArcOriginalNewHS::~ArcOriginalNewHS() {
  if (ofile_arcboundsOriginalDat != NULL) fclose(ofile_arcboundsOriginalDat);
  delete[] emitted_flat_evap;
}

void ArcOriginalNewHS::print_par() const {
  printf( " - beta_tip:            %g \n",          this->beta_tip         );
  printf( " - beta_f:              %g \n",          this->beta_f           );
  printf( " - j_melt:              %g [A/cm^2]\n",  this->j_melt           );
  double j_max = 4.7133e9 * SQU(12.0) * exp(-62.338/12.0);
  if (j_melt > j_max) {
    printf( "   (NOTE: j_melt > j_max = %g [A/cm^2] -- No melting possible!)\n", j_max);
  }
  printf( " - do_erosion:          %c \n", this->do_erosion ? 'y' : 'n'    );
  printf( " - Remission:           %g [dz]\n",      this->Remission        );
  printf( " - Remission_theor:     %g [dz]\n",      this->Remission_theor  );
  printf( " - evap_flat_center:    %c \n", this->evap_flat_center ? 'y' : 'n');
  printf( " - Rborder:             %g [dz]\n",      this->Rborder          );
  printf( " - SEY:                 %g \n",          this->SEY              );
  printf( " - r_Cu_e:              %g \n",          this->r_Cu_e           );
  printf( " - r_Cu_e_flat:         %g \n",          this->r_Cu_e_flat      );
  printf( " - evap_flat_center:    %c \n", this->evap_flat_center ? 'y' : 'n');
  printf( " - alpha_flat:          %g \n",          this->alpha_flat       );
  printf( " - fracInjectStep:      %c \n", this->fracInjectStep ? 'y' : 'n');
  printf( " - heatspike_threshold: %g [particles/cm^2/s]\n",  this->heatspike_threshold);
  printf( "   = %g [particles/Ldb^2/injection_timestep]\n",
	  this->heatspike_threshold*(Ndb*n2inj_step*Omega_pe) /
	  (4.1938226e7*n_ref*sqrt(T_ref)) );
  printf( "   = %g [particles/inner cell/injection_timestep]\n",
	  this->heatspike_threshold*(Ndb*n2inj_step*Omega_pe) /
	  (4.1938226e7*n_ref*sqrt(T_ref)) * PI*SQU(dr) );
  printf( "   = %g [particles/2nd cell/injection_timestep]\n",
	  this->heatspike_threshold*(Ndb*n2inj_step*Omega_pe) /
	  (4.1938226e7*n_ref*sqrt(T_ref)) * 3*PI*SQU(dr) );
  printf( " - heatspike_yieldRatio: %g [outgoing/incomming]\n", this->heatspike_yield_ratio);
  printf( " - heatspike_model:     %c \n",          this->heatspike_model  );
  printf( " - file_timestep:       %u \n",          this->file_timestep    );
}

void ArcOriginalNewHS::init(unsigned int nr, double zmin, double zmax, double rmax) {
  ArcBounds::init(nr,zmin,zmax,rmax);
  
  v_inj_e = v_te*0.01;
  v_inj_i = vt_ions[2];

  sput_cathode_current = new unsigned int[nr];
  for (unsigned int i = 0; i < nr; i++) {
    sput_cathode_current[i] = 0;
  }

  //Initialization
  this->has_melted          = false;
  this->erosion             = 0;
  this->emitted_tip_melting = 0;
  
  this->emitted_tip_evap    = 0;
  this->emitted_flat_evap   = new int[nr];
  for (unsigned int i = 0; i < nr; i++) {
    emitted_flat_evap[i] = 0;
  }

  this->emitted_tip_output         = 0;
  this->emitted_flat_output        = 0;
  this->emitted_SEY_output         = 0;
  this->emitted_evap_output        = 0;
  this->emitted_sputter_cat_output = 0;
  this->emitted_sputter_ano_output = 0;
  this->emitted_heatspike_output   = 0;

  this->heatspike_sigma    = 0;
  this->heatspike_incident = 0;

}
void ArcOriginalNewHS::re_init(unsigned int nr, double zmin, double zmax, double rmax) {
  //ArcBounds::re_init(nr, zmin, zmax, rmax);
  ArcOriginalNewHS::init(nr, zmin, zmax, rmax); //ArcBounds::re_init just calls init()
}

void ArcOriginalNewHS::remove_i(Particle pa[], size_t &np, unsigned int sort) {
  if(sort != 1){
    if (np != 0) {
      printf("Error detected in ArcOriginal::remove_i(): %zu particles of sort=%u\n", np,sort);
      exit(1);
    }
    return;
  }

  size_t n_lost = 0;

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

      Sput newSput = calc_sput(pa[n], cs_ions[sort], sput_cathode_current);
      if ( newSput.Y != 0 ) {
	sput_cathode.push_back(newSput);
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
  
}
void ArcOriginalNewHS::remove_n(Particle pa[], size_t &np){
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

      Sput newSput = calc_sput(pa[n], cs_ions[NSpecies-1], sput_cathode_current);
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

Sput ArcOriginalNewHS::calc_sput(const Particle& pa, const double cs, unsigned int* current_enhancedY) {
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

void ArcOriginalNewHS::inject_e(Particle pa[], size_t &np, double const Ez[]) {
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

  if (Eloc > 0.) {  
    // Protect against numerical fluctuations
    // (limit of FN formula validity)
    if ( Eloc > 12. ) Eloc = 12.;
    
    j = 4.7133e9 * SQU(Eloc) * exp(-62.338/Eloc); // in A/cm^2
    j *= PI*(Ndb*e2inj_step*Omega_pe*SQU(Remission_theor*dz))/(6.7193e-12*n_ref*sqrt(T_ref)); //dimless, 2D: rem^2
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
      if ( Rborder <= jj ) 
	jFN[jj] *= (Ndb*e2inj_step*Omega_pe*PI*SQU(dz)*(2*jj+1))/(6.7193e-12*n_ref*sqrt(T_ref));
      else if ( (Rborder > jj) && (Rborder < jj+1) )
	jFN[jj] *= Ndb*e2inj_step*Omega_pe*PI*SQU(dz)*( SQU(jj+1)-SQU(Rborder) )/(6.7193e-12*n_ref*sqrt(T_ref));
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
	printf("Error in ArcOriginalNewHS::inject_e(): Particle array overflow (SEY)\n");
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
  emitted_SEY_output += totY;
  sput_cathode_SEY.clear();

  // FN INJECTION FROM FIELD EMITTER
  // FE: inject with flat distribution over the emission radius
  //Inject 'tmp' num particles
  size_t tmp = size_t(j);
  j -= tmp;
  if ( RAND <= j ) tmp++;
  
  if ( ( np + tmp ) >= NPART) {
    printf("Error in ArcOriginalNewHS::inject_e(): Particle array overflow (FN at field emitter)\n");
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
  
  injected_e[0]       += tmp;  
  emitted_tip_evap    += tmp;
  emitted_tip_melting += tmp;
  emitted_tip_output  += tmp;
  
  // FN INJECTION OUTSIDE THE FIELD EMITTER
  for (unsigned int jj=0; jj<nr; jj++ ) {
    if ( jFN[jj] > 0.0 ) {
      //Inject 'tmp' num. particles
      tmp = (size_t) (jFN[jj]); 
      jFN[jj] -= tmp;
      if ( RAND <= jFN[jj] ) tmp++;
      
      if ( ( np + tmp ) >= NPART) {
	printf("Error in ArcOriginalNewHS::inject_e(): Particle array overflow (FN outside field emitter)\n");
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
	  if ( Rborder <= jj ) {
	    //pa[np+k].p.r = jj + RAND;
	    pa[np+k].p.r = sqrt( (2*jj+1)*RAND + SQU(jj) );
	  }
	  else if ( (Rborder > jj) && (Rborder < jj+1) ) {
	    //pa[np+k].p.r = Rborder + (jj + 1 - Rborder) * RAND;
	    pa[np+k].p.r = sqrt( ( SQU(jj+1)-SQU(Rborder) )*RAND + SQU(Rborder) );
	  }
	  else {
	    printf("ERROR in ArcOriginalNewHS::inject_e():\n");
	    printf("unexpected non-zero jFN - jj=%u, Rborder=%g\n", jj, Rborder);
	    printf("nsteps=%d jFN[jj]=%g\n, k=%zu, pa[np+k].p.r=%g\n\n", nsteps, jFN[jj]+tmp, k, pa[np+k].p.r); //nsteps = global variable!	    
	  }

	  pa[np+k].p.z = zmin;
	  
	  //Fractional timestep push (euler method)
	  double Rp = RAND; //1-R; how (fractionally) far into the timestep
	  // are we at z=0?
	  pa[np+k].p.r  += Rp*e2inj_step*pa[np+k].p.vr;
	  pa[np+k].p.z  += Rp*e2inj_step*pa[np+k].p.vz;
	  //No magnetic field, E = Ez on surface
	  pa[np+k].p.vz -= (Rp*e2inj_step-0.5)*2*Ez[(jj+1)*NZ];
	  
	  //Reflect on axis
	  if ( pa[np+k].p.r < 0 ) pa[np+k].p.r = -pa[np+k].p.r;
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

	  if ( Rborder <= jj ) {
	    //pa[np+k].p.r = jj + RAND + pa[np+k].p.vr;
	    pa[np+k].p.r = sqrt( (2*jj+1)*RAND + SQU(jj) ) + pa[np+k].p.vr;
	  }
	  else if ( (Rborder > jj) && (Rborder < jj+1) ) {
	    //pa[np+k].p.r = Rborder + (jj + 1 - Rborder) * RAND + pa[np+k].p.vr;
	    pa[np+k].p.r = sqrt( ( SQU(jj+1)-SQU(Rborder) )*RAND + SQU(Rborder) ) + pa[np+k].p.vr;
	  }
	  else {
	    printf("ERROR in ArcOriginalNewHS::inject_e():\n");
	    printf("unexpected non-zero jFN - jj=%u, Rborder=%g\n", jj, Rborder);
	    printf("nsteps=%d jFN[jj]=%g\n, k=%zu, pa[np+k].p.r=%g\n\n", nsteps, jFN[jj]+tmp, k, pa[np+k].p.r); //nsteps = global variable!	    
	  }

	  if ( pa[np+k].p.r < 0 ) pa[np+k].p.r = 1.e-20;
	  else if (pa[np+k].p.r > nr) pa[np+k].p.r = 2*nr - pa[np+k].p.r;
	}
	pa[np+k].p.m = 1;
	current_e[0] += 2*( int(pa[np+k].p.r) ) + 1;
	current_cathode[ int(pa[np+k].p.r) ] += 1;
      }
      np += tmp;

      injected_e[0]         += tmp;
      emitted_flat_evap[jj] += tmp;
      emitted_flat_output   += tmp;
    }
  } 
}

void ArcOriginalNewHS::inject_n(Particle pa[], size_t &np, double const Ez[]) {

  // Cathode, neutral evaporation on tip
  double tmp = 0;
  if (evap_flat_center) {
    int sum_flat = 0;
    for (unsigned int i = 0; i < nr; i++) {
      sum_flat += emitted_flat_evap[i];
      emitted_flat_evap[i] = 0;
    }
    tmp = r_Cu_e * emitted_tip_evap + r_Cu_e_flat * sum_flat;
    emitted_tip_evap = 0;
  }
  else {
    tmp = r_Cu_e * emitted_tip_evap;
    emitted_tip_evap = 0;
  }
  size_t n2inject_evap = (size_t) tmp;
  tmp -= n2inject_evap;
  if ( RAND <= tmp ) n2inject_evap++;

  double r1, r2;
  if ( ( np + n2inject_evap ) >= NPART) {
    printf("Error in ArcOriginalNewHS::inject_n(): Particle array overflow (tip evaporation)\n");
    exit(1);
  }
  for ( size_t k=0; k<n2inject_evap; k++ ) {
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
  np                   += n2inject_evap;
  erosion              += n2inject_evap;
  emitted_evap_output  += n2inject_evap;
  injected_n[0]        += n2inject_evap;
  
  // Cathode, neutral evaporation on flat surface
  if (not evap_flat_center) {
    for (unsigned int i = 0; i < nr; i++) {
      tmp = r_Cu_e * emitted_flat_evap[i];
      emitted_flat_evap[i] = 0;
  
      n2inject_evap = (size_t)tmp;
      tmp -= n2inject_evap;
      if ( RAND <= tmp ) n2inject_evap++;
    
      if ( ( np + n2inject_evap ) >= NPART) {
	printf("Error in ArcOriginalNewHS::inject_n(): Particle array overflow (flat evaporation)\n");
	exit(1);
      }
      for ( size_t k=0; k<n2inject_evap; k++ ) {
	// Velocity
	do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
	r2 = RAND * TWOPI;
	pa[np+k].p.vr = r1*cos(r2)*v_inj_i; // do not suppress
	pa[np+k].p.vt = r1*sin(r2)*v_inj_i; // do not suppress
	
	do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
	pa[np+k].p.vz = r1*v_inj_i;
	
	// Position
	pa[np+k].p.z = zmin + pa[np+k].p.vz;

	if ( Rborder <= i ) {
	  //pa[np+k].p.r = i + RAND + pa[np+k].p.vr;
	  pa[np+k].p.r = sqrt( (2*i+1)*RAND + SQU(i) ) + pa[np+k].p.vr;
	}
	else if ( (Rborder > i) && (Rborder < i+1) ) {
	  //pa[np+k].p.r = Rborder + (i + 1 - Rborder) * RAND + pa[np+k].p.vr;
	  pa[np+k].p.r = sqrt( ( SQU(i+1)-SQU(Rborder) )*RAND + SQU(Rborder) ) + pa[np+k].p.vr;
	}
	
	if ( pa[np+k].p.r < 0 ) pa[np+k].p.r = 1.e-20;
	else if (pa[np+k].p.r > nr) pa[np+k].p.r = 2*nr - pa[np+k].p.r;
	
	pa[np+k].p.m = 1; //SN;
	
      }
      np                  += n2inject_evap;
      erosion             += n2inject_evap;
      emitted_evap_output += n2inject_evap;
      injected_n[0]       += n2inject_evap;
    }
  }
  
  // Cathode, Yamamura&Tawara-sputtering
  size_t foo = np;
  inject_sput(pa, np, sput_cathode, true);
  emitted_sputter_cat_output += (np-foo);

  // Cathode, heat-spike sputtering

  //Threshold for heatspike sputtering [superparticles /  lambda_db^2 / t_inj_neutrals
  double threshold = heatspike_threshold*(Ndb*n2inj_step*Omega_pe) / 
    (4.1938226e7*n_ref*sqrt(T_ref));
  
  if (heatspike_model == 'a' or heatspike_model == 'b') {
    unsigned int num_incident        = 0; //Number of superparticles in cells > threshold
    unsigned int num_incident_inside = 0; //Number of superparticles in cells since last threshold crossing
    unsigned int num_incident2       = 0; //Used together with num_incident_inside
    unsigned int sigma               = 0; //Sigma of radial gaussian distribution
                                          // to inject superparticles in case of heatspike

    for (unsigned int i=0; i<nr; i++ ) {
      unsigned int cell_incident = sput_cathode_current[i];
      // Rescale current with area of cell / area Ldb^2
      double current = cell_incident / ( PI * (2*i + 1) * SQU(dr) );
      
      num_incident_inside += cell_incident;
      
      if ( current >= threshold ) {
	num_incident += cell_incident;
	
	num_incident2 += num_incident_inside;
	num_incident_inside = 0;
	
	sigma = i+1;
      }
      
      sput_cathode_current[i] = 0;
    }

    //store for output and sputtering
    heatspike_sigma = sigma;
    if (heatspike_model == 'b') heatspike_incident = num_incident2;
    else /* a */                heatspike_incident = num_incident;
  }
  else if (heatspike_model == 'c') {
    unsigned int num_incident = 0;
    unsigned int num_incident_inside = 0;
    unsigned int sigma = 0;

    for (unsigned int i = 0; i<nr; i++) {
      unsigned int cell_incident = sput_cathode_current[i];
      num_incident += cell_incident;
      double current = num_incident / ( PI*SQU((i+1)*dr) );
      
      if ( current >= threshold ) {
	num_incident_inside = num_incident;
	sigma = i+1;
      }
      sput_cathode_current[i] = 0;
    }
    
    heatspike_sigma = sigma;
    heatspike_incident = num_incident_inside;
  }
  else {
    printf("Error: No valid heatspike_model specified ('%c')\n", heatspike_model);
    exit(1);
  }

  //If we are above the threshold, sputter.
  if ( heatspike_incident > 0 ) {
    // Enhanced yield from MD
    // Generate Gaussian distribution for r-coord's 
    double num_sputtered = heatspike_incident * heatspike_yield_ratio;
    //Handle floating point part
    unsigned int num_sputtered_i = (unsigned int) ( num_sputtered );
    if (RAND <= num_sputtered - num_sputtered_i) {
      num_sputtered_i += 1;
    }
    double r = 0; //Temp variable
    for (unsigned int j=0; j<num_sputtered_i; j++ ) {
      do { r = heatspike_sigma * sqrt(-2.*log(RAND+1.e-20)); } while ( r >= nr );
      Sput foo;
      foo.r = r;
      foo.Y = 1;
      sput_cathode.push_back(foo);
    }
    inject_sput(pa, np, sput_cathode, true);
    emitted_heatspike_output += num_sputtered_i;
  }
  
  // Anode, sputtering only
  foo = np;
  inject_sput(pa, np, sput_anode, false);
  emitted_sputter_ano_output += (np-foo);

}
void ArcOriginalNewHS::inject_sput(Particle pa[], size_t &np, vector<Sput> &sput, bool isCathode) {
  
  double r1, r2;

  size_t tot = 0;
  for ( size_t k=0; k<sput.size(); k++ ) {
    for ( int i=0; i<sput[k].Y; i++ ) {
      if ( ( np + tot ) >= NPART) {
	printf("Error in ArcOriginalNewHS::inject_sput(): Particle array overflow\n");
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


void ArcOriginalNewHS::timestep(unsigned int nstep, bool isOutputTimestep) {
  if ( !has_melted ) {
    if (do_erosion) {
      beta_tip -= erosion * n_ref / ( Ndb * pow(Remission_theor*dz,3)) * 3.824131e-24; 
      if (beta_tip < beta_f) { beta_tip = beta_f; }
    }
    this->erosion = 0;
  }

  // Melting of the emitter due to Ohmic heating by FE current
  if (not has_melted && emitted_tip_melting > 
      ( j_melt*Ndb*Omega_pe*SQU(Remission_theor*dz)*PI/6.7193e-12/n_ref/sqrt(T_ref) ) ) {
    beta_tip = beta_f;
    has_melted = true;
  }
  emitted_tip_melting = 0;
    
  //Write files and reset counting arrays
  ArcBounds::timestep(nstep, isOutputTimestep);
}

void ArcOriginalNewHS::writeFile_arcboundsOriginalDat(unsigned int nstep) {
  if (ofile_arcboundsOriginalDat == NULL) {
    ofile_arcboundsOriginalDat = fopen("arcbounds_original.dat", "w");
    fprintf(ofile_arcboundsOriginalDat, "## timestep beta_tip has_melted[y/n] emitted_tip[numSP/dtOut] E_tip_loc[GV/m] emitted_flat[numSP/dtOut] emit_SEY[numSP/dtOut] emit_evap[numSP/dtOut] emit_sput_cat[numSP/dtOut] emit_sput_ano[numSP/dtOut] emit_htspk[numSP/dtOut] htspk_sig[dz] htspk_inc[numSP/dt]\n");
    fflush(ofile_arcboundsOriginalDat);
  }

  if ( nstep % file_timestep == 0) {
    if (nstep % n2inj_step == 0) {
      fprintf(ofile_arcboundsOriginalDat, "%11d %f %c %d %f %d %d %d %d %d %d %u %u\n", 
	      nstep, beta_tip, (has_melted? 'y' : 'n'), emitted_tip_output, E_tip_loc, 
	      emitted_flat_output, emitted_SEY_output, emitted_evap_output,
	      emitted_sputter_cat_output, emitted_sputter_ano_output, emitted_heatspike_output,
	      heatspike_sigma, heatspike_incident);
    }
    else {
      fprintf(ofile_arcboundsOriginalDat, "%11d %f %c %d %f %d %d -1 -1 -1 -1 -1 -1\n",
	      nstep, beta_tip, (has_melted? 'y' : 'n'),  emitted_tip_output, E_tip_loc,
	      emitted_flat_output, emitted_SEY_output);
    }
    fflush(ofile_arcboundsOriginalDat);
    
    emitted_tip_output         = 0;
    emitted_flat_output        = 0;
    emitted_SEY_output         = 0;
    emitted_evap_output        = 0;
    emitted_sputter_cat_output = 0;
    emitted_sputter_ano_output = 0;
    emitted_heatspike_output   = 0;
  }
}
