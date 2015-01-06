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

  arcbounds_flexFN2.cpp:
  Particle boundary conditions,
  Fowler-Nordheim where beta and alpha are functions of the r-coordinate

***********************************************************************/

#include "arcbound_flexFN.h"
#include "mydef.h"

#include "random.h"

#define   XTRN extern
#include "var.h"
#include "arrays1.h"
#undef XTRN

#include <cmath>

using namespace std;

// ******** Implementation of FlexFN ******************
void FlexFN::calcFN_current(double const Ez[], double* alpha, double* beta,
		    double* currentFN) {
  for (size_t i = 0; i <= nr; i++) { //Loop over meshpoints
    double field = Ez[i*NZ];
    double deltaI_FN = 0.0;
    if (field < 0) {
      // Fowler-Nordheim with Wang-Loew approximation
      // W-L: v(y) ~ 0.956 - 1.062*(3.7947e-5)^2*Eloc/(4.5*4.5)
      // work fct=4.5eV, j in A/cm^2, E in V/m
      
      // rescale the field to GV/m, multiply with beta(r_i)
      field = - 2.69036254e-10*dz/SQU(Omega_pe)*sqrt(T_ref*n_ref)*field*beta[i];
      // Protect against numerical fluctuations
      // (this caps j at the field where tunneling probability = 1 (check!)
      if ( field > 12. ) field = 12.;      
      deltaI_FN = 4.7133e9 * SQU(field) * exp(-62.338/field); // in A/cm^2
      
      //Rescale to units (#Superparticles / omega_pe^-1) / lambda_Db^2
      deltaI_FN *= Ndb/(6.7192539e-12*n_ref*sqrt(T_ref));
      //Area factor
      deltaI_FN *= alpha[i];
      //Find the current (#Superparticles) through each annular disk
      // between (j - 0.5)*DZ and (j+0.5)*DZ AND inside the domain (0, nr*dz)
      if (i == 0) 
	deltaI_FN *= (PI*SQU(0.5*dz))*(e2inj_step*Omega_pe);
      else if (i > 0 && i < nr)
	deltaI_FN *= (TWOPI*i*SQU(dz))*(e2inj_step*Omega_pe);
      else if (i == nr)
	deltaI_FN *= (PI*(nr-0.25)*SQU(dz))*(e2inj_step*Omega_pe);
    }
    currentFN[i] += deltaI_FN;    
  } // END loop over meshpoints
}

void FlexFN::injectFN(Particle pa[], size_t &np, double* currentFN, double const Ez[]) {

  double v_inj_e = v_te*0.01;

  for (size_t i = 0; i <= nr; i++) { //Loop over meshpoints
    
    //Number to inject
    size_t Ninj = size_t (currentFN[i]);
    if ( RAND <= currentFN[i]-Ninj ) Ninj++;
    //Safety check
    if (Ninj + np > NPART) {
      printf("Error in FlexFN::injectFN(): Particle array overflow (FN at field emitter)\n");
      exit(1);
    }
    
    //Inject!
    double r1; //Temp variable
    for (size_t k = 0; k < Ninj; k++) {
      
      //Velocity
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
      double r2 = RAND * TWOPI;
      pa[np+k].p.vr = r1*cos(r2)*v_inj_e;
      pa[np+k].p.vt = r1*sin(r2)*v_inj_e;
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
      pa[np+k].p.vz = r1*v_inj_e;
      
      //Position
      // uniform random position on annulus
      double a = i > 0  ? (i-0.5) : 0.0;
      double b = i < nr ? (i+0.5) : nr;
      pa[np+k].p.r = sqrt((b*b-a*a)*RAND+a*a);
      pa[np+k].p.z = zmin;
      
      //Fractional timestep push (euler method)
      double Rp = RAND; //1-R; how (fractionally) far into the timestep
                        // are we at z=0?
      pa[np+k].p.r  += Rp*e2inj_step*pa[np+k].p.vr;
      pa[np+k].p.z  += Rp*e2inj_step*pa[np+k].p.vz;
      //Assume no magnetic field, E = Ez on surface
      pa[np+k].p.vz -= (Rp*e2inj_step-0.5)*2*Ez[i*NZ];
      
      if ( pa[np+k].p.r < 0 ) pa[np+k].p.r = -pa[np+k].p.r; //Reflect on axis
      else if (pa[np+k].p.r > nr ) pa[np+k].p.r = 2*nr - pa[np+k].p.r;
      
      pa[np+k].p.m = 1;
      
      current_e[0] += 2*( int(pa[np+k].p.r) ) + 1;
      current_cathode[ int(pa[np+k].p.r) ] += 1;
    }

    np += Ninj;
    injected_e[0] += Ninj;
    
  } // END loop over meshpoints
}

// ******** Implementation of FlexFN_ring ******************
FlexFN_ring::FlexFN_ring(std::vector<char*>& options) {
  if (options.size() != 5) {
    cout << "Error in FlexFN_ring(): Expected 5 options, got "
	 << options.size() << endl;
    exit(1);
  }

  sscanf(options[0], "%*[^:]%*[:] %lg", &(this->alpha_ring));
  sscanf(options[1], "%*[^:]%*[:] %lg", &(this->beta_ring));
  sscanf(options[2], "%*[^:]%*[:] %u",  &(this->idx1_ring));
  sscanf(options[3], "%*[^:]%*[:] %u",  &(this->idx2_ring));
  sscanf(options[4], "%*[^:]%*[:] %u",  &(this->file_timestep));
}
void FlexFN_ring::print_par() const {
  printf( " - alpha_ring           %g \n", alpha_ring);
  printf( " - beta_ring            %g \n", beta_ring);
  printf( " - idx1_ring            %u \n", idx1_ring);
  printf( " - idx2_ring            %u \n", idx2_ring);
  printf( " - file_timestep        %u \n", file_timestep );
}
FlexFN_ring::~FlexFN_ring() {
  delete[] FN_alpha;
  delete[] FN_beta;
  delete[] FN_current;
}

void FlexFN_ring::init(unsigned int nr, double zmin, double zmax, double rmax) {
  ArcBounds::init(nr,zmin,zmax,rmax);
  FN_alpha   = new double[nr+1];
  FN_beta    = new double[nr+1];
  FN_current = new double[nr+1];
  initAlphaBeta();
}
void FlexFN_ring::re_init(unsigned int nr, double zmin, double zmax, double rmax) {
  FlexFN_ring::init(nr, zmin, zmax, rmax);
}

void FlexFN_ring::initAlphaBeta() {
  for (unsigned int i = 0; i < nr+1; i++) {
    if (i >= idx1_ring && i <= idx2_ring) {
      FN_alpha[i] = alpha_ring;
      FN_beta[i]  = beta_ring;
    }
    else {
      FN_alpha[i] = 0.0;
      FN_beta[i]  = 0.0;
    }
  }
}

void FlexFN_ring::inject_e(Particle pa[], size_t &np, double const Ez[]) {
  for (unsigned int i = 0; i < nr+1; i++) FN_current[i] = 0;
  calcFN_current(Ez,FN_alpha,FN_beta,FN_current);
  injectFN(pa, np, FN_current, Ez);
}

// ******** Implementation of FlexFN_twoComp ******************
FlexFN_twoComp::FlexFN_twoComp(std::vector<char*>& options) {
  if (options.size() != 6) {
    cout << "Error in FlexFN_twoComp(): Expected 6 options, got "
	 << options.size() << endl;
    exit(1);
  }

  sscanf(options[0], "%*[^:]%*[:] %lg", &(this->alpha1));
  sscanf(options[1], "%*[^:]%*[:] %lg", &(this->alpha2));
  sscanf(options[2], "%*[^:]%*[:] %lg", &(this->beta1));
  sscanf(options[3], "%*[^:]%*[:] %lg", &(this->beta2));
  sscanf(options[4], "%*[^:]%*[:] %u",  &(this->idx1));
  sscanf(options[5], "%*[^:]%*[:] %u",  &(this->file_timestep));
}
void FlexFN_twoComp::print_par() const {
  printf( " - alpha1               %g \n", alpha1);
  printf( " - alpha2               %g \n", alpha2);
  printf( " - beta1                %g \n", beta1);
  printf( " - beta2                %g \n", beta2);
  printf( " - idx1                 %u \n", idx1);
  printf( " - file_timestep        %u \n", file_timestep );
}

void FlexFN_twoComp::init(unsigned int nr, double zmin, double zmax, double rmax) {
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

}
void FlexFN_twoComp::re_init(unsigned int nr, double zmin, double zmax, double rmax) {
  FlexFN_twoComp::init(nr, zmin, zmax, rmax);
}
FlexFN_twoComp::~FlexFN_twoComp() {
  for (size_t i = 0; i < 2; i++) {
    delete[] FN_alpha[i];
    delete[] FN_beta[i];
  }
  delete[] FN_alpha;
  delete[] FN_beta;

  delete[] FN_current;
}


void FlexFN_twoComp::inject_e(Particle pa[], size_t &np, double const Ez[]) {
  for (unsigned int i = 0; i < nr+1; i++) FN_current[i] = 0;
  
  for (unsigned int i = 0; i < 2; i++) {
    calcFN_current(Ez,FN_alpha[i],FN_beta[i],FN_current);
  }
  
  injectFN(pa, np, FN_current, Ez);
}

