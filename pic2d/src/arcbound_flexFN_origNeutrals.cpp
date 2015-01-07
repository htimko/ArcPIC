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

  arcbounds_flexFN_origNeutrals.cpp:  
  Particle boundary conditions,
  Fowler-Nordheim where beta and alpha are functions of the r-coordinate,
  with sputtering and evaporation similar to Helga's thesis.

***********************************************************************/

#include "arcbound_flexFN_origNeutrals.h"
#include "mydef.h"

#include "random.h"

#define   XTRN extern
#include "var.h"
#include "arrays1.h"
#undef XTRN

#include <cmath>

using namespace std;

// ******** Implementation of FlexFN_twoComp_origNeutrals ******************
FlexFN_twoComp_origNeutrals::FlexFN_twoComp_origNeutrals(std::vector<char*>& options) {
  if (options.size() != 9) {
    cout << "Error in FlexFN_twoComp_origNeutrals(): Expected 9 options, got "
	 << options.size() << endl;
    exit(1);
  }

  sscanf(options[0], "%*[^:]%*[:] %lg", &(this->alpha1));
  sscanf(options[1], "%*[^:]%*[:] %lg", &(this->alpha2));
  sscanf(options[2], "%*[^:]%*[:] %lg", &(this->beta1));
  sscanf(options[3], "%*[^:]%*[:] %lg", &(this->beta2));
  sscanf(options[4], "%*[^:]%*[:] %u",  &(this->idx1));
  sscanf(options[5], "%*[^:]%*[:] %lg", &(this->SEY));
  sscanf(options[6], "%*[^:]%*[:] %lg", &(this->r_Cu_e));
  
  char foo;
  sscanf(options[7],  "%*[^:]%*[:] %c",  &(foo));
  if      (foo == 'y') this->doHeatspike = true;
  else if (foo == 'n') this->doHeatspike = false;
  else {
    printf("Error in FlexFN_twoComp_origNeutrals::FlexFN_twoComp_origNeutrals(): doHeatspike has to be either 'y' or 'n' \n");
    exit(1);
  }

  sscanf(options[8], "%*[^:]%*[:] %u",  &(this->file_timestep));
}
void FlexFN_twoComp_origNeutrals::print_par() const {
  printf( " - alpha1               %g \n", alpha1);
  printf( " - alpha2               %g \n", alpha2);
  printf( " - beta1                %g \n", beta1);
  printf( " - beta2                %g \n", beta2);
  printf( " - idx1                 %u \n", idx1);
  printf( " - SEY:                 %g \n", SEY);
  printf( " - r_Cu_e:              %g \n", r_Cu_e);
  printf( " - doHeatspike          %c \n", doHeatspike ? 'y' : 'n');
  printf( " - file_timestep        %u \n", file_timestep );
}

void FlexFN_twoComp_origNeutrals::init(unsigned int nr, double zmin, double zmax, double rmax) {
  ArcBounds::init(nr,zmin,zmax,rmax);
  
  FN_alpha    = new double*[2];
  FN_alpha[0] = new double[nr+1];
  FN_alpha[1] = new double[nr+1];

  FN_beta     = new double*[2];
  FN_beta[0]  = new double[nr+1];
  FN_beta[1]  = new double[nr+1];

  for (unsigned int i = 0; i < nr+1; i++) {
    if (i <= idx1) {
      FN_alpha[0][i] = alpha1;
      FN_beta[0][i]  = beta1;
    }
    else {
      FN_alpha[0][i] = 0.0;
      FN_beta[0][i]  = 0.0;
    }
    FN_alpha[1][i] = alpha2;
    FN_beta[1][i]  = beta2;
  }

  FN_current = new double[nr+1];
  
  FN_current_sum = new double[nr+1];
  for (unsigned int i = 0; i < nr+1; i++) FN_current_sum[i] = 0.0;


}
void FlexFN_twoComp_origNeutrals::re_init(unsigned int nr, double zmin, double zmax, double rmax) {
  FlexFN_twoComp_origNeutrals::init(nr, zmin, zmax, rmax);
}
FlexFN_twoComp_origNeutrals::~FlexFN_twoComp_origNeutrals() {
  for (size_t i = 0; i < 2; i++) {
    delete[] FN_alpha[i];
    delete[] FN_beta[i];
  }
  delete[] FN_alpha;
  delete[] FN_beta;

  delete[] FN_current;
  
  delete[] FN_current_sum;
}


void FlexFN_twoComp_origNeutrals::inject_e(Particle pa[], size_t &np, double const Ez[]) {
  //Fowler-nordheim two-component injection
  for (unsigned int i = 0; i < nr+1; i++) FN_current[i] = 0;
  for (unsigned int i = 0; i < 2; i++) {
    calcFN_current(Ez,FN_alpha[i],FN_beta[i],FN_current);
  }
  injectFN(pa, np, FN_current, Ez);
  for (unsigned int i = 0; i < nr+1; i++) FN_current_sum[i] += FN_current[i];

  //SEY at cathode
  double r1, r2; // temp variables from Gaussian RNG
  double v_inj_e = v_te*0.01; //electron injection velocity
  // SEY: inject with incident r-coordinates, empty SEY variables!
  size_t totY = 0;
  for (size_t k=0; k < sput_cathode_SEY.size(); k++ ) {
    for (int i=0; i < sput_cathode_SEY[k].Y; i++ ) {
      if ( ( np + totY ) >= NPART ) {
	printf("Error in FlexFN_twoComp_origNeutrals::inject_e(): Particle array overflow (SEY)\n");
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
      else if ( pa[np+totY].p.r > nr) pa[np+totY].p.r = 2*nr - pa[np+totY].p.r;

      pa[np+totY].p.m = 1;
      
      current_e[0] += 2*( int(pa[np+totY].p.r) ) + 1;

      totY++;
    }
  }
  np += totY;
  injected_e[0] += totY;
  sput_cathode_SEY.clear();
}

void FlexFN_twoComp_origNeutrals::inject_n(Particle pa[], size_t &np,
					   double const Ez[]) {

  double v_inj_i = vt_ions[2];

  // Cathode, neutral evaporation (per-cell)
  double r1, r2;
  for (size_t i = 0; i <= nr; i++) { //Loop over meshpoints/injection cells
    //Calculate number to inject
    double tmp = r_Cu_e * FN_current_sum[i];
    size_t n2inject_evap = size_t ( tmp );
    tmp -= n2inject_evap;
    if ( RAND <= tmp ) n2inject_evap++;
    FN_current_sum[i] = 0.0;
    
    if ( ( np + n2inject_evap ) >= NPART) {
      printf("Error in FlexFN_twoComp_origNeutrals::inject_n(): Particle array overflow ( evaporation in cell %zu )\n", i);
      exit(1);
    }
    
    //Inject the particles: Uniform distribution over the cell
    for ( size_t k=0; k<n2inject_evap; k++ ) {
      // Velocity:
      // Gaussian scheme
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
      r2 = RAND * TWOPI;
      pa[np+k].p.vr = r1*cos(r2)*v_inj_i; // do not suppress
      pa[np+k].p.vt = r1*sin(r2)*v_inj_i; // do not suppress
      // Gaussian scheme
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
      pa[np+k].p.vz = r1*v_inj_i;
      
      // uniform random position on annulus
      double a = i > 0  ? (i-0.5) : 0.0;
      double b = i < nr ? (i+0.5) : nr;
      pa[np+k].p.r = sqrt((b*b-a*a)*RAND+a*a);
      pa[np+k].p.z = zmin;
      
      //Fractional timestep push (euler method, no acceleration)
      double Rp = RAND; //1-R; how (fractionally) far into the timestep
                        // are we at z=0?
      pa[np+k].p.r  += Rp*n2inj_step*pa[np+k].p.vr;
      pa[np+k].p.z  += Rp*n2inj_step*pa[np+k].p.vz;
      
      if ( pa[np+k].p.r < 0 ) pa[np+k].p.r = -pa[np+k].p.r; //Reflect on axis
      else if ( pa[np+k].p.r > nr) pa[np+k].p.r = 2*nr - pa[np+k].p.r;

      pa[np+k].p.m = 1; //SN;    
    }
    np += n2inject_evap;
    injected_n[0] += n2inject_evap;
  }
  
  // Cathode, sputtering
  inject_sput(pa, np, sput_cathode, true, v_inj_i);
  // Anode, sputtering only
  inject_sput(pa, np, sput_anode, false, v_inj_i);
}
void FlexFN_twoComp_origNeutrals::inject_sput(Particle pa[], size_t &np, std::vector<Sput> &sput, bool isCathode, double v_inj) {
  
  double r1, r2;

  size_t m = sput.size();
  size_t tot = 0;
  for ( size_t k=0; k<m; k++ ) {
    for ( int i=0; i<sput[k].Y; i++ ) {
      if ( ( np + tot ) >= NPART) {
	printf("Error in FlexFN_twoComp_origNeutrals::inject_n(): Particle array overflow (sputtering)\n");
	exit(1);
      }
      
      // Gaussian scheme      
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
      r2 = RAND * TWOPI;
      pa[np+tot].p.vr = r1*cos(r2)*v_inj;
      pa[np+tot].p.vt = r1*sin(r2)*v_inj; 
      
      // Gaussian scheme
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
      pa[np+tot].p.vz = isCathode ? r1*v_inj : -r1*v_inj;
      
      pa[np+tot].p.z = (isCathode ? zmin : zmax) + pa[np+tot].p.vz;
      pa[np+tot].p.r = sput[k].r + pa[np+tot].p.vr; 
      
      if ( pa[np+tot].p.r < 0 ) pa[np+tot].p.r = 1.e-20;
      else if ( pa[np+tot].p.r > nr) pa[np+tot].p.r = 2*nr - pa[np+tot].p.r;

      pa[np+tot].p.m = 1; //SN;
      
      // Sum up total yield for outputting
      tot++;
    }
  }
  np += tot;
  injected_n[isCathode ? 0 : 1] += tot;
  sput.clear();
}
void FlexFN_twoComp_origNeutrals::remove_i(Particle pa[], size_t &np, unsigned int sort) {
  if(sort != 1){
    if (np != 0) {
      printf("Error detected in FlexFN_twoComp_origNeutrals::remove_i(): %zu particles of sort=%u\n", np,sort);
      exit(1);
    }
    return;
  }
  
  size_t n_lost = 0;

  // check threshold of sputtering yield, only for ions, only at the cathode
  // threshold = 1e7 A/cm^2
  double threshold = 1.e7*(Ndb*dt_ion*Omega_pe) / 
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
  
  if ( check_enh == false  or doHeatspike == false) {
    // Yamamura-Tawara fitting for Cu -> Cu
    for (size_t j=0; j<sput_cathode_temp.size(); j++ ) {
      sput_cathode.push_back(sput_cathode_temp[j]);
    }
  }
  else if ( check_enh == true ) {
    // Enhanced yield from MD
    // Generate Gaussian distribution for r-coord's 
#warning "Strange model - always generate 1000 Cu's?"
    for (int j=0; j<1000; j++ ) { // 1000/SN
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

  if (check_enh == true and doHeatspike == false) {
    printf("Skipped heatspike sputtering as it is switched off, yield=%d\n", Ycat_sum);
    fflush(stdout);
  }
}
void FlexFN_twoComp_origNeutrals::remove_n(Particle pa[], size_t &np){
  size_t n_lost = 0; 
  
  for (size_t n=0; n<np; n++ ) {
    //"Infinity"
    if ( pa[n].p.r >= rmax ) {
      removedNeutrals.push_back(pa[n]);
      removed_n[2]++;
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

Sput FlexFN_twoComp_origNeutrals::calc_sput(const Particle& pa,
					     const double cs,
					     double* current_enhancedY) {
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
