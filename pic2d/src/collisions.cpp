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

  collisions.cpp:
  2D (r,z) collision routines after Konstantin Matyash
    x -> z 
    y -> r
    z -> t

***********************************************************************/

#include <cmath>
#include <cstdlib>
#include <cstdio>

#include <vector>

#include <omp.h>

#include  "pic.h"
#include  "dim.h"
#include  "mydef.h"
#include  "random.h"

#include "colls.h"

#define   XTRN extern
#include  "var.h"

/**********************************************************************

 e-e or i-i Coulomb collision routine (for the same type of particles)
 kind = 0 e-e
 otherwise i-i

 Input:
 - ordcount        Ordering array, ordcount[i] = number of particles
                     in cell i=ir*NZ+iz
 - nr              Number of cells in r-direction
 - nz              Number of cells in z-direction
 - NZ              Number of grids in z-direction = nz+1
 - Mpa_over_me     Mass of colliding particles/electron mass
 - kind            Type of particle: 0 for electrons, != 0 for ions
 - ncoll           How many timesteps between each call to collisions
 In/Out:
 - Particle pa[]   Array of particles to collide, sorted by cell
 Out:
 - momcheck        Diagnostics, checking that momentum is conserved
 - engcheck        Diagnostics, checking that energy is conserved

***********************************************************************/

void coll_el_knm_2D( Particle  pa[], size_t ordcount[],
		     int nr,int nz,int NZ, 
		     double Mpa_over_me, int kind,
		     Vec3d *momcheck, double *engcheck, int ncoll) {
  
  Vec3d	v_rel, v_rel_DELTA;
  
  static std::vector<size_t> inds2coll; //Not nice for parallelization
                                        // (neither over spieces- or cell)

  static const double Amplcoulomb =  1.;  // Amplification of coulomb collisions, only for testing purposes
  static const double LanLog = 13.;       // Coulomb Log 
  //Constant factor in coulomb collisions
  double Acoll = Amplcoulomb * ( !kind ? 1.0 : dt_ion*dt_ion*dt_ion ) * LanLog * SQU(SQU(Omega_pe)) * ncoll / 
    ( TWOPI*PI * SQU(SQU(dz))*SQU(dz) * SQU(Mpa_over_me) * SQU(Ndb) * N_sp );
  
  size_t Next = 0; // Index in pa[] of first particle in cell ir*NZ+iz (the way this is done today messes up paralellization)
  for (int ir=0; ir<nr;  ir++) {
    for (int iz=0; iz<nz; iz++) { 	        

      //Number of particles left to collide in cell (ir,iz)
      size_t N2coll = ordcount[ir*NZ+iz];
      //Resize the indices
      inds2coll.resize(N2coll);

      //Constant factor in columb collisions for this cell
      double Acoll_cell = Acoll*ordcount[ir*NZ+iz]/(2.0*ir+1.0);

      if ( N2coll>1 ) {
	//Pick two random particle indices (j,k) from the particles that has not yet collided
	for (size_t loopind  = 0; loopind < N2coll; loopind++ ) inds2coll[loopind]=loopind;
	
	size_t j(0), k(0);
	while ( N2coll>0 ) {
	  if ( N2coll>1 ) {
	    size_t jind = (size_t)(N2coll*RAND);
	    if (jind == N2coll) jind--;   
	    j = inds2coll[jind];
	    for (size_t loopind = jind + 1; loopind < N2coll; loopind++) inds2coll[loopind-1] = inds2coll[loopind];
	    N2coll--;
	    
	    jind = (size_t)(N2coll*RAND);
	    if (jind == N2coll) jind--;     
	    k = inds2coll[jind];
	    for (size_t loopind = jind + 1; loopind < N2coll; loopind++) inds2coll[loopind-1] = inds2coll[loopind];
	    N2coll--;
	  }
	  else {
	    //Last particle in a cell with odd number of particles: collide it with "j" from last time
	    k = inds2coll[N2coll-1];
	    N2coll--;
	  }
	  
	  // BEGIN Just for diagnostic purposes - remove when done
	  (*momcheck).z += pa[j + Next].p.vz + pa[k + Next].p.vz;
	  (*momcheck).r += pa[j + Next].p.vr + pa[k + Next].p.vr;
	  (*momcheck).t += pa[j + Next].p.vt + pa[k + Next].p.vt;
	  *engcheck += SQU(pa[j + Next].p.vz) + SQU(pa[j + Next].p.vr) + SQU(pa[j + Next].p.vt) + 
	    SQU(pa[k + Next].p.vz) + SQU(pa[k + Next].p.vr) + SQU(pa[k + Next].p.vt);
	  // END

	  // relative velocity
	  v_rel.z = pa[j + Next].p.vz - pa[k + Next].p.vz;
	  v_rel.r = pa[j + Next].p.vr - pa[k + Next].p.vr;
	  v_rel.t = pa[j + Next].p.vt - pa[k + Next].p.vt;
	  // W = ||v_rel|| (and inverse)
	  double W  = sqrt( SQU(v_rel.z) + SQU(v_rel.r) + SQU(v_rel.t));
	  //printf("\tW=%e\n",W);
	  if (W < 1.e-10) {
	    //if u->0, the relative change might be big, but it's a big change on top of nothing => SKIP
	    continue; 
	  }
	  double IW = 1./W;
	  
	  //MC to find the scattering angle, through transformation of delta
	  // which is istributed with P(delta) = N(0, <delta^2>)
	  double deltaVar2 = Acoll_cell*IW*IW*IW; // Variance <delta^2>
	  double delta = GausRandom(0,sqrt(deltaVar2));
	  double sinTheta = 2*delta/(1+SQU(delta));
	  double oneMinusCosTheta = sinTheta*delta;
	  //Scattering polar angle
	  double Phi = RAND*TWOPI;
	  double cosPhi = cos(Phi);
	  double sinPhi = sin(Phi);


	  //v_rel projected into Z-R plane (and inverse)
	  double vpZR  = sqrt( SQU(v_rel.z) + SQU(v_rel.r) );
	  if (vpZR > 1e-10) {
	    //Normal case, need rotation of coordinate system etc.
	    double ivpZR = 1./vpZR;
	    v_rel_DELTA.z = v_rel.z*v_rel.t*ivpZR*cosPhi*sinTheta - v_rel.z*oneMinusCosTheta - W*v_rel.r*ivpZR*sinPhi*sinTheta;
	    v_rel_DELTA.r = ivpZR*sinTheta*(v_rel.t*v_rel.r*cosPhi + v_rel.z*W*sinPhi) - v_rel.r*oneMinusCosTheta;
	    v_rel_DELTA.t = -vpZR*cosPhi*sinTheta - v_rel.t*oneMinusCosTheta;
	  }
	  else {
	    //v_rel in t-direction only => coordinate systems are identical => simpler equations
	    v_rel_DELTA.z = W*sinTheta*cosPhi;
	    v_rel_DELTA.r = W*sinTheta*sinPhi;
	    v_rel_DELTA.t = -W*oneMinusCosTheta;
	  }

	  //Update the particle velocities
	  pa[j + Next].p.vz += v_rel_DELTA.z/2.;    
	  pa[j + Next].p.vr += v_rel_DELTA.r/2.;
	  pa[j + Next].p.vt += v_rel_DELTA.t/2.;
	  
	  pa[k + Next].p.vz -= v_rel_DELTA.z/2.;
	  pa[k + Next].p.vr -= v_rel_DELTA.r/2.;
	  pa[k + Next].p.vt -= v_rel_DELTA.t/2.;
	  
	  // BEGIN Just for diagnostic purposes - remove it when U done
	  (*momcheck).z -= pa[j + Next].p.vz + pa[k + Next].p.vz;
	  (*momcheck).r -= pa[j + Next].p.vr + pa[k + Next].p.vr;
	  (*momcheck).t -= pa[j + Next].p.vt + pa[k + Next].p.vt;
	  *engcheck -= SQU(pa[j + Next].p.vz) + SQU(pa[j + Next].p.vr) + SQU(pa[j + Next].p.vt) + SQU(pa[k + Next].p.vz) + SQU(pa[k + Next].p.vr) + SQU(pa[k + Next].p.vt);
	  // END (of diagnostics )

	}//END do-loop over particles
      }//END if (N2coll > 1)
      Next+=ordcount[ir*NZ+iz];
      //END loop over cells
    }
  }
}


/**********************************************************************

 Ion - neutrals elastic charge exchange, no Super neutrals

***********************************************************************/

void coll_ion_neutral_noSP_2D( Particle  neutrals[], size_t ordcount_ntrls[], double M_n,
			       Particle  ions[], size_t ordcount_ion[], double M_i,
			       int nr, int nz, int NZ, Reaction React,
			       Vec3d *momcheck, double *engcheck   ) {

  double W2_0, W_0, IW_0, v_proj_zr, ivp, cos_theta, sin_theta, cos_beta, sin_beta;
  double cos_phi, sin_phi, psi;

  double vz, vr, vt; 
  Vec3d	v_rel;                                               

  static std::vector<size_t> indxs; //Not nice for parallelization
                                    // (neither over spieces- or cell)

  //Scattering fit
  double Emin  = React.Emin;
  //double Emax  = React.Emax;
  double Estep = React.Estep;
  int Epoints = React.N;
  int Eind;
  double S_i = 0.;
      
  size_t Next_ion = 0;
  size_t Next_n  = 0;
  
  double dti  = 1./dt_ion;  // Scaling factor 4 ion velocity

  for (int jr=0; jr<nr; jr++) {
    for (int jz=0; jz<nz; jz++) {
      size_t n_n  = ordcount_ntrls[jr*NZ+jz];
      size_t n_ion = ordcount_ion[jr*NZ+jz];

      if ( n_n && n_ion ) {
	if ( n_ion < n_n ) {
	  
	  size_t N2coll = n_n;  // only for SP = 1 !!   
	  indxs.resize(n_ion);
	  
	  for (size_t i_ion = 0; i_ion < n_ion; i_ion++ ) {
	    size_t k = (size_t)(N2coll*RAND);
	    if (k == N2coll) k--;                            		         
	    indxs[i_ion] = k;  
	  }
	  
	  for (size_t i_ion = 0; i_ion < n_ion; i_ion++) {
	    size_t i_n = indxs[i_ion];
	    
	    // relative velocity
	    v_rel.z = ions[i_ion + Next_ion].p.vz - neutrals[i_n + Next_n].p.vz;
	    v_rel.r = ions[i_ion + Next_ion].p.vr - neutrals[i_n + Next_n].p.vr;       
	    v_rel.t = ions[i_ion + Next_ion].p.vt - neutrals[i_n + Next_n].p.vt;
	    
	    // before scattering
	    W2_0 = SQU(v_rel.z) + SQU(v_rel.r) + SQU(v_rel.t);               
	    W_0  = sqrt(W2_0);
	    
	    // Linear fit for Cross-section
	    Eind = (int)((W2_0 - Emin)/Estep);
	    if ( Eind < 0)
	      S_i = React.CS[0];
	    else if ( Eind >= Epoints)
	      S_i = React.CS[Epoints];
	    else
	      S_i = React.CS[Eind] + (React.CS[Eind +1] - React.CS[Eind])*((W2_0 - Emin)/Estep - Eind);
	    S_i *= Ampl;

	    // 1-exp(-n_n*W_0*Si*dt_coll)
	    if ( RAND < n_n*W_0*S_i*dti/(2*jr+1) ) { // correct Si with 2D volume factor (2j+1) here   
	      *engcheck += 0.5*M_i*(SQU(ions[i_ion + Next_ion].p.vz) + SQU(ions[i_ion + Next_ion].p.vr) + SQU(ions[i_ion + Next_ion].p.vt));
	      *engcheck += 0.5*M_n*(SQU(neutrals[i_n + Next_n].p.vz) + SQU(neutrals[i_n + Next_n].p.vr) + SQU(neutrals[i_n + Next_n].p.vt));
	      
	      (*momcheck).z += M_i*ions[i_ion + Next_ion].p.vz + M_n*neutrals[i_n + Next_n].p.vz;
	      (*momcheck).r += M_i*ions[i_ion + Next_ion].p.vr + M_n*neutrals[i_n + Next_n].p.vr;
	      (*momcheck).t += M_i*ions[i_ion + Next_ion].p.vt + M_n*neutrals[i_n + Next_n].p.vt;
	      
	      W_0  = MAX( W_0, 1.e-10);   
	      IW_0 = 1./W_0;                          
	      
	      v_proj_zr = sqrt( SQU(v_rel.z) + SQU(v_rel.r) ); 
	      v_proj_zr = MAX( v_proj_zr, 1.e-10);              
	      ivp       = 1./v_proj_zr;
	      
	      cos_theta = v_rel.t * IW_0;     
	      sin_theta = v_proj_zr*IW_0;
	      cos_beta = v_rel.z*ivp;    
	      sin_beta = v_rel.r*ivp;
	      
	      if (RAND < 0.5) {
#warning What is this?!? (found multiple places)
		// just 4 test!! Use = (2.+E-2.*pow((1+E), RAND))/E instead
		cos_phi = 2.*RAND -1. ;
	      }
	      else {
		cos_phi = -1.;
	      }
	      
	      sin_phi = sqrt(1. - SQU(cos_phi));
	      psi  = TWOPI*RAND;
	      
	      vt = W_0*cos_phi;     
	      vz = W_0*sin_phi*cos(psi);
	      vr = W_0*sin_phi*sin(psi);
	      
	      v_rel.z -= vz*cos_beta*cos_theta - vr*sin_beta + vt*cos_beta*sin_theta;
	      v_rel.r -= vz*sin_beta*cos_theta + vr*cos_beta + vt*sin_beta*sin_theta;  
	      v_rel.t -= vt*cos_theta - vz*sin_theta;
	      
	      ions[i_ion + Next_ion].p.vz -= M_n/(M_n + M_i)*v_rel.z;
	      ions[i_ion + Next_ion].p.vr -= M_n/(M_n + M_i)*v_rel.r;
	      ions[i_ion + Next_ion].p.vt -= M_n/(M_n + M_i)*v_rel.t;
	      
	      neutrals[i_n + Next_n].p.vz    += M_i/(M_n + M_i)*v_rel.z; 
	      neutrals[i_n + Next_n].p.vr    += M_i/(M_n + M_i)*v_rel.r;     
	      neutrals[i_n + Next_n].p.vt    += M_i/(M_n + M_i)*v_rel.t;    
	      
	      *engcheck -= 0.5*M_i*(SQU(ions[i_ion + Next_ion].p.vz) + SQU(ions[i_ion + Next_ion].p.vr) + SQU(ions[i_ion + Next_ion].p.vt));
	      *engcheck -= 0.5*M_n*(SQU(neutrals[i_n + Next_n].p.vz) + SQU(neutrals[i_n + Next_n].p.vr) + SQU(neutrals[i_n + Next_n].p.vt));
	      
	      (*momcheck).z -= M_i*ions[i_ion + Next_ion].p.vz + M_n*neutrals[i_n + Next_n].p.vz;
	      (*momcheck).r -= M_i*ions[i_ion + Next_ion].p.vr + M_n*neutrals[i_n + Next_n].p.vr;
	      (*momcheck).t -= M_i*ions[i_ion + Next_ion].p.vt + M_n*neutrals[i_n + Next_n].p.vt;
	      
	    }   // if (W2_0  RAND < n_n*W_0*Si)
	    
	  }   // for ( i_ion = 0; i_ion < n_ion; i_ion++)
	}   // if (n_ion < n_n)

	/* Case when n_n <= n_ion*/       
	else {
	  /********************************************************
	       Let's make array of random indexes for ions here
	  ********************************************************/
	  
	  size_t N2coll = n_ion;
	  indxs.resize(n_n);

	  for (size_t i_n = 0; i_n < n_n; i_n++) {
	    size_t k = (size_t)(N2coll*RAND);
	    if (k == N2coll) k--;
	    indxs[i_n] = k;
	  }
           
	  for (size_t i_n = 0; i_n < n_n; i_n++) {
	    size_t i_ion = indxs[i_n];
	    
	    // relative velocity
	    v_rel.z= ions[i_ion + Next_ion].p.vz - neutrals[i_n + Next_n].p.vz;  
	    v_rel.r= ions[i_ion + Next_ion].p.vr - neutrals[i_n + Next_n].p.vr;          
	    v_rel.t= ions[i_ion + Next_ion].p.vt - neutrals[i_n + Next_n].p.vt;
	    
	    // before scattering
	    W2_0 = SQU(v_rel.z) + SQU(v_rel.r) + SQU(v_rel.t);  
	    W_0  = sqrt(W2_0);
	    
	    // Linear fit for Cross-section
	    Eind = (int)((W2_0 - Emin)/Estep);
	    if ( Eind < 0)
	      S_i = React.CS[0];
	    else if ( Eind >= Epoints)
	      S_i = React.CS[Epoints];
	    else
	      S_i = React.CS[Eind] + (React.CS[Eind +1] - React.CS[Eind])*((W2_0 - Emin)/Estep - Eind);
	    S_i *= Ampl;
	    
	    
	    // 1-exp(-n_n*W_0*Si*dt_coll)
	    if ( RAND < n_ion*W_0*S_i*dti/(2*jr+1) ) { // correct Si with 2D volume factor (2j+1) here     
	      *engcheck += 0.5*M_i*(SQU(ions[i_ion + Next_ion].p.vz) + SQU(ions[i_ion + Next_ion].p.vr) + SQU(ions[i_ion + Next_ion].p.vt));
	      *engcheck +=  0.5*M_n*(SQU(neutrals[i_n + Next_n].p.vz) + SQU(neutrals[i_n + Next_n].p.vr) + SQU(neutrals[i_n + Next_n].p.vt));
	      
	      (*momcheck).z += M_i*ions[i_ion + Next_ion].p.vz + M_n*neutrals[i_n + Next_n].p.vz;
	      (*momcheck).r += M_i*ions[i_ion + Next_ion].p.vr + M_n*neutrals[i_n + Next_n].p.vr;
	      (*momcheck).t += M_i*ions[i_ion + Next_ion].p.vt + M_n*neutrals[i_n + Next_n].p.vt;
	      
	      W_0  = MAX( W_0, 1.e-10);       
	      IW_0 = 1./W_0;                          
	      
	      v_proj_zr = sqrt( SQU(v_rel.z) + SQU(v_rel.r) );   
	      v_proj_zr = MAX( v_proj_zr, 1.e-10);              
	      ivp       = 1./v_proj_zr;
	      
	      cos_theta = v_rel.t * IW_0;   
	      sin_theta = v_proj_zr*IW_0;
	      cos_beta = v_rel.z*ivp;  
	      sin_beta = v_rel.r*ivp;
	      
	      
	      if (RAND < 0.5) {
		// just 4 test!! Use = (2.+E-2.*pow((1+E), RAND))/E instead
		cos_phi = 2.*RAND -1. ; 
	      }
	      else {
		cos_phi = -1.;
	      }
	      
	      sin_phi = sqrt(1. - SQU(cos_phi));
	      psi  = TWOPI*RAND;
	      
	      vt = W_0*cos_phi;          
	      vz = W_0*sin_phi*cos(psi);
	      vr = W_0*sin_phi*sin(psi);
	      
	      v_rel.z -= vz*cos_beta*cos_theta - vr*sin_beta + vt*cos_beta*sin_theta;  
	      v_rel.r -= vz*sin_beta*cos_theta + vr*cos_beta + vt*sin_beta*sin_theta;   
	      v_rel.t -= vt*cos_theta - vz*sin_theta;
	      
	      ions[i_ion + Next_ion].p.vz -= M_n/(M_n + M_i)*v_rel.z;
	      ions[i_ion + Next_ion].p.vr -= M_n/(M_n + M_i)*v_rel.r;
	      ions[i_ion + Next_ion].p.vt -= M_n/(M_n + M_i)*v_rel.t;
	      
	      neutrals[i_n + Next_n].p.vz += M_i/(M_n + M_i)*v_rel.z;   
	      neutrals[i_n + Next_n].p.vr += M_i/(M_n + M_i)*v_rel.r;   
	      neutrals[i_n + Next_n].p.vt += M_i/(M_n + M_i)*v_rel.t;    
	      
	      
	      *engcheck -= 0.5*M_i*(SQU(ions[i_ion + Next_ion].p.vz) + SQU(ions[i_ion + Next_ion].p.vr) + SQU(ions[i_ion + Next_ion].p.vt));
	      *engcheck -= 0.5*M_n*(SQU(neutrals[i_n + Next_n].p.vz) + SQU(neutrals[i_n + Next_n].p.vr) + SQU(neutrals[i_n + Next_n].p.vt));
	      
	      (*momcheck).z -= M_i*ions[i_ion + Next_ion].p.vz + M_n*neutrals[i_n + Next_n].p.vz;
	      (*momcheck).r -= M_i*ions[i_ion + Next_ion].p.vr + M_n*neutrals[i_n + Next_n].p.vr;
	      (*momcheck).t -= M_i*ions[i_ion + Next_ion].p.vt + M_n*neutrals[i_n + Next_n].p.vt;
	      
	    }   // if (W2_0  RAND < n_n*W_0*Si)            
	    
	  }   // for ( i_n = 0; i_n < n_n; i_n++)
	  
	  
	}   //  if n_ion < n_n else
      }   //if ( n_n && n_ion )
      
      Next_ion += n_ion;
      Next_n  += n_n;
      
    }
  }  //for (j=0; j < Ngrid;  j++)

}




/**********************************************************************

 Neutral- neutrals ellastic, no Super neutrals

***********************************************************************/

void coll_n_n_2D( Particle neutrals[], size_t ordcount_ntrls[],                              
		  int nr, int nz, int NZ, Reaction React,
		  Vec3d *momcheck, double *engcheck ) {
  
  //Setup domain decomposition
  // 1. count total number of particles && make cumulative ordcount
  static size_t* ordcount_cum = NULL;
  if (ordcount_cum == NULL) ordcount_cum = new size_t[nr*nz+1];   
  ordcount_cum[0] = 0;
  
  for (int cidx=0; cidx<nr*nz; cidx++) {
    //Linear continious cell indexing
    //Double loop easier, but this is a good testbed
    // for paralellizable single loop
    int jr = cidx/nz;
    int jz = cidx-jr*nz;
    //printf("%i/%i , %i/%i\n", jr, nr, jz, nz);
    //np += ordcount_ntrls[jr*NZ+jz];
    ordcount_cum[cidx+1] = ordcount_cum[cidx] + ordcount_ntrls[jr*NZ+jz];
  }
  size_t np = ordcount_cum[nr*nz];

  // 2. Distribute the cells
  //  This array is used for load balancing,
  //  contains the linear continous cell index indicating
  //  where the CPU should start processing.
  static size_t* splitIdx = NULL;
  if (splitIdx == NULL) splitIdx = new size_t[numParaThreads+1];
  //for (int i = 0; i < numParaThreads; i++) splitIdx[i] = 0;
  splitIdx[numParaThreads] = nr*nz;

  if (np <= COLL_N_N_2D_OMP_MINPARTICLES or numParaThreads == 1) {
    //Too few particles: don't split
    splitIdx[0] = 0;
    for (int i = 1; i < numParaThreads; i++) splitIdx[i] = nr*nz;
  }
  else {
    int particles_per_CPU = np / numParaThreads;
    
    splitIdx[0] = 0;
    int orderIdx = 1;

    int np_now = 0;
    for (int cidx=0; cidx<nr*nz; cidx++) {
      //Linear continious cell indexing
      int jr = cidx/nz;
      int jz = cidx-jr*nz;
      np_now += ordcount_ntrls[jr*NZ+jz];
      if (np_now > particles_per_CPU) {
	splitIdx[orderIdx] = cidx;
	orderIdx += 1;
	np_now = 0;
      }
    }
    //Make sure the whole splitIdx array is initialized
    for (; orderIdx < numParaThreads; orderIdx++) splitIdx[orderIdx] = nr*nz;
  }

  //Initialize sorting arrays
  static std::vector<size_t>* sub_2coll_all = NULL;
  static std::vector<size_t>* indxs_all     = NULL;
  if (sub_2coll_all == NULL) sub_2coll_all = new std::vector<size_t> [numParaThreads];
  if (indxs_all     == NULL) indxs_all     = new std::vector<size_t> [numParaThreads];
     
  //Here comes the paralellization stuff!
#pragma omp parallel if(splitIdx[1] != ((size_t)nr)*((size_t)nz)) 
  {
    
    int CPUIDX = omp_get_thread_num();
    
    double W2_0, W_0, IW_0, v_proj_zr, ivp, cos_theta, sin_theta, cos_beta, sin_beta;
    double cos_phi, sin_phi, psi;
    
    double vz, vr, vt; 
    
    Vec3d v_rel;                                               
    
    double one_or_half;  // half if without pair            
    
    double Emin  = React.Emin;
    //double Emax  = React.Emax;
    double Estep = React.Estep;
    int  Epoints = React.N;
    int  Eind;
    double S_i = 0.;
    
    size_t Next_n  = ordcount_cum[splitIdx[CPUIDX]];
    
    double dti  = 1./dt_ion;  // Scaling factor for the ion velocity
    
    std::vector<size_t>& sub_2coll = sub_2coll_all[CPUIDX];
    std::vector<size_t>& indxs     = indxs_all    [CPUIDX];
    
    for (size_t cidx=splitIdx[CPUIDX]; cidx<splitIdx[CPUIDX+1]; cidx++) {      
      //Linear continious cell indexing
      int jr = cidx/nz;
      int jz = cidx-jr*nz;
      
      size_t n_n  = ordcount_ntrls[jr*NZ+jz];
      sub_2coll.resize(n_n);
      indxs.resize(n_n);
      
      // At least one pair
      if ( n_n > 1 ) {
	
	size_t N2coll = n_n;     // only for SP = 1 !!
	for (size_t i=0; i<n_n;  i++ ) {
	  sub_2coll[i] = 1;      // same here	    
	}
	
	/********************************************************
	 Let's make array of random indexes for neutrals here
	********************************************************/
	
	for (size_t i=0; i<n_n; i++ ) { 
	  size_t k = (size_t)(N2coll*Random(CPUIDX));
	  if (k == N2coll) k--;
	  
	  size_t i_n    = 0;             
	  size_t i_next = sub_2coll[i_n];
	  
	  while ( i_next <= k ) {
	    i_next += sub_2coll[++i_n];
	  }
	  
	  --sub_2coll[i_n];
	  N2coll--;
	  
	  indxs[i] = i_n;
	}
	
	
	for (size_t i = 0; i < n_n; i += 2) {
	  size_t i_n1 = indxs[i];
	  size_t i_n2 = 0;
	  if ( i+1 < n_n ){
	    i_n2 = indxs[i+1];
	    one_or_half = 1.;
	  }
	  else {
	    i_n2 = indxs[1];
	    // odd number - last one gets collided every second time (on average)
	    one_or_half = 0.5; 
	  }  
	  
	  // relative velocity
	  v_rel.z = neutrals[i_n1 + Next_n].p.vz - neutrals[i_n2 + Next_n].p.vz; 
	  v_rel.r = neutrals[i_n1 + Next_n].p.vr - neutrals[i_n2 + Next_n].p.vr;
	  v_rel.t = neutrals[i_n1 + Next_n].p.vt - neutrals[i_n2 + Next_n].p.vt;
	  
	  // before scattering
	  W2_0 = SQU(v_rel.z) + SQU(v_rel.r) + SQU(v_rel.t);               
	  W_0  = sqrt(W2_0);
	  
	  // Linear fit for Cross-section
	  Eind = (int)((W2_0 - Emin)/Estep);
	  if ( Eind < 0)
	    S_i = React.CS[0];
	  else if ( Eind >= Epoints)
	    S_i = React.CS[Epoints];
	  else
	    S_i = React.CS[Eind] + (React.CS[Eind +1] - React.CS[Eind])*((W2_0 - Emin)/Estep - Eind);
	  S_i *= Ampl*one_or_half;
	  
	  // 1-exp(-n_n*W_0*Si*dt_coll)
	  if ( Random(CPUIDX) < n_n*W_0*S_i*dti/(2*jr+1) ) { // correct Si with 2D volume factor (2j+1) here
	    
	    *engcheck += SQU(neutrals[i_n1 + Next_n].p.vz) + SQU(neutrals[i_n1 + Next_n].p.vr) + SQU(neutrals[i_n1 + Next_n].p.vt);
	    *engcheck += SQU(neutrals[i_n2 + Next_n].p.vz) + SQU(neutrals[i_n2 + Next_n].p.vr) + SQU(neutrals[i_n2 + Next_n].p.vt);
	    
	    (*momcheck).z += neutrals[i_n1 + Next_n].p.vz + neutrals[i_n2 + Next_n].p.vz;
	    (*momcheck).r += neutrals[i_n1 + Next_n].p.vr + neutrals[i_n2 + Next_n].p.vr;
	    (*momcheck).t += neutrals[i_n1 + Next_n].p.vt + neutrals[i_n2 + Next_n].p.vt;
	    
	    W_0  = MAX( W_0, 1.e-10);       // e-10 should be small enough
	    IW_0 = 1./W_0;
	    
	    v_proj_zr = sqrt( SQU(v_rel.z) + SQU(v_rel.r) );
	    v_proj_zr = MAX( v_proj_zr, 1.e-10);              
	    ivp       = 1./v_proj_zr;
	    
	    cos_theta = v_rel.t * IW_0;              
	    sin_theta = v_proj_zr*IW_0;
	    cos_beta = v_rel.z*ivp;               
	    sin_beta = v_rel.r*ivp;
	    
	    cos_phi = 2.*Random(CPUIDX) -1. ; 
	    sin_phi = sqrt(1. - SQU(cos_phi));
	    psi  = TWOPI*Random(CPUIDX);
	    
	    vt = W_0*cos_phi;             
	    vz = W_0*sin_phi*cos(psi);
	    vr = W_0*sin_phi*sin(psi);
	    
	    v_rel.z -= vz*cos_beta*cos_theta - vr*sin_beta + vt*cos_beta*sin_theta;    
	    v_rel.r -= vz*sin_beta*cos_theta + vr*cos_beta + vt*sin_beta*sin_theta;   
	    v_rel.t -= vt*cos_theta - vz*sin_theta;
	    
	    neutrals[i_n1 + Next_n].p.vz -= 0.5*v_rel.z;
	    neutrals[i_n1 + Next_n].p.vr -= 0.5*v_rel.r;
	    neutrals[i_n1 + Next_n].p.vt -= 0.5*v_rel.t;
	    
	    neutrals[i_n2 + Next_n].p.vz += 0.5*v_rel.z;  
	    neutrals[i_n2 + Next_n].p.vr += 0.5*v_rel.r;    
	    neutrals[i_n2 + Next_n].p.vt += 0.5*v_rel.t;    
	    
	    *engcheck -= SQU(neutrals[i_n1 + Next_n].p.vz) + SQU(neutrals[i_n1 + Next_n].p.vr) + SQU(neutrals[i_n1 + Next_n].p.vt);
	    *engcheck -= SQU(neutrals[i_n2 + Next_n].p.vz) + SQU(neutrals[i_n2 + Next_n].p.vr) + SQU(neutrals[i_n2 + Next_n].p.vt);
	    
	    (*momcheck).z -= neutrals[i_n1 + Next_n].p.vz + neutrals[i_n2 + Next_n].p.vz;
	    (*momcheck).r -= neutrals[i_n1 + Next_n].p.vr + neutrals[i_n2 + Next_n].p.vr;
	    (*momcheck).t -= neutrals[i_n1 + Next_n].p.vt + neutrals[i_n2 + Next_n].p.vt;
	    
	    
	  }   //  if (W2_0  RAND < n_n*W_0*Si)
	  
	}   //for ( i_ion = 0; i_ion < n_ion; i_ion++)
      }   //if ( n_n > 1 )
      
      Next_n  += n_n;
      
    } //END for (cidx -> jr, jz)
  } //END omp parallel
}




/**********************************************************************

 Simplified version of inelastic collisions of electrons with whatever - 
 we just let the electron loose some energy and then rotate it's velocity
 reaction products are not taken into account

***********************************************************************/

void coll_el_all_fake_2D( Particle molecules[], size_t ordcount_m[], double M_m,   // molecules
			  Particle electrons[], size_t ordcount_el[],              // electrons
			  int nr, int nz, int NZ, Reaction React) {
  
  const double m_e = 1.;
  
  double W2, W, IW, W2_0, W_0, IW_0, v_proj_zr, ivp, cos_theta, sin_theta, cos_beta, sin_beta;
  //double W0;
  double cos_phi, sin_phi, psi;
  
  double vz, vr, vt;
  Vec3d v_rel, delta_v_rel, target_v;       
  
  double E_th  = React.Eth;
  double Emin  = React.Emin;
  //double Emax  = React.Emax;
  double Estep = React.Estep;
  int Epoints = React.N;
  int Eind;
  double S_i = 0.;
  
  size_t Next_el = 0;
  size_t Next_m  = 0;
  
  double dti  = 1./dt_ion;    // Scaling factor 4 ion velocity
  
  for (int jr=0; jr<nr; jr++) {
    for (int jz=0; jz<nz; jz++) {
      size_t n_m  = ordcount_m[jr*NZ+jz];
      size_t n_el = ordcount_el[jr*NZ+jz];
      
      if ( n_m && n_el ) {
	size_t nm_real  = 0;
	
	for (size_t i=0; i<n_m; i++) {
	  nm_real += molecules[i + Next_m].p.m;            
	}
	
	//size_t N2coll = nm_real;              
	
	for (size_t i_el=0; i_el<n_el; i_el++) {
	  size_t i_m = (size_t)(n_m*RAND);
	  if (i_m == n_m) i_m--;              
	  
	  // relative velocity
	  v_rel.z = electrons[i_el + Next_el].p.vz - molecules[i_m + Next_m].p.vz*dti; 
	  v_rel.r = electrons[i_el + Next_el].p.vr - molecules[i_m + Next_m].p.vr*dti;
	  v_rel.t = electrons[i_el + Next_el].p.vt - molecules[i_m + Next_m].p.vt*dti;
	  
	  // before scattering
	  W2_0 = SQU(v_rel.z) + SQU(v_rel.r) + SQU(v_rel.t);                
	  W_0  = sqrt(W2_0);
	  
	  // Linear fit for Cross-section
	  if ( W2_0 >= E_th ) {
	    Eind = (int)((W2_0 - Emin)/Estep);
	    if ( Eind < 0)
	      S_i = React.CS[0];
	    else if ( Eind >= Epoints)
	      S_i = React.CS[Epoints];
	    else
	      S_i = React.CS[Eind] + (React.CS[Eind +1] - React.CS[Eind])*((W2_0 - Emin)/Estep - Eind);
	  }
	  S_i *= Ampl;
	  
	  // 1-exp(-n_m*W_0*Si*dt_coll)
	  if (W2_0 >= E_th && RAND < nm_real*W_0*S_i/(2*jr+1) ) { // correct Si with 2D volume factor (2j+1) here 
	    
	    /*
	     *engcheck += 0.5*m_e*(SQU(electrons[i_el + Next_el].p.vz) + SQU(electrons[i_el + Next_el].p.vr) + SQU(electrons[i_el + Next_el].p.vt));
	     *engcheck += 0.5*M_m*(SQU(molecules[i_m + Next_m].p.vz*dti) + SQU(molecules[i_m + Next_m].p.vr*dti) + SQU(molecules[i_m + Next_m].p.vt*dti)) - 0.5*m_e*M_m/(m_e+M_m)*E_th ;
	     
	     (*momcheck).z += m_e*electrons[i_el + Next_el].p.vz + M_m*molecules[i_m + Next_m].p.vz*dti;
	     (*momcheck).r += m_e*electrons[i_el + Next_el].p.vr + M_m*molecules[i_m + Next_m].p.vr*dti;
	     (*momcheck).t += m_e*electrons[i_el + Next_el].p.vt + M_m*molecules[i_m + Next_m].p.vt*dti;
	    */
	    
	    W_0  = MAX( W_0, 1.e-14);   
	    IW_0 = 1./W_0;
	    
	    W2 = W2_0 - E_th;
	    W  = sqrt(W2);    
	    
	    // After energy loss
	    delta_v_rel.z = v_rel.z*(W*IW_0 -1.);
	    delta_v_rel.r = v_rel.r*(W*IW_0 -1.);     
	    delta_v_rel.t = v_rel.t*(W*IW_0 -1.);
	    
	    v_rel.z += delta_v_rel.z;
	    v_rel.r += delta_v_rel.r;
	    v_rel.t += delta_v_rel.t;
	    
	    W  = MAX( W, 1.e-10);    
	    IW = 1./W;
	    
	    
	    v_proj_zr = sqrt( SQU(v_rel.z) + SQU(v_rel.r) ); 
	    v_proj_zr = MAX( v_proj_zr, 1.e-10);              
	    ivp       = 1./v_proj_zr;
	    
	    cos_theta = v_rel.t * IW;   
	    sin_theta = v_proj_zr*IW;
	    cos_beta = v_rel.z*ivp;                
	    sin_beta = v_rel.r*ivp;
	    
	    
	    /* Sampling scattering angles in C.M. system */
	    // just 4 test!! Use = (2.+E-2.*pow((1+E), RAND))/E instead
	    cos_phi = 2.*RAND -1. ; 
	    sin_phi = sqrt(1. - SQU(cos_phi));
	    psi  = TWOPI*RAND;
	    
	    vt = W*cos_phi;           
	    vz = W*sin_phi*cos(psi);
	    vr = W*sin_phi*sin(psi);
	    
	    v_rel.z -= vz*cos_beta*cos_theta - vr*sin_beta + vt*cos_beta*sin_theta; 
	    v_rel.r -= vz*sin_beta*cos_theta + vr*cos_beta + vt*sin_beta*sin_theta;   
	    v_rel.t -= vt*cos_theta - vz*sin_theta;
	    
	    v_rel.z -= delta_v_rel.z;
	    v_rel.r -= delta_v_rel.r;
	    v_rel.t -= delta_v_rel.t;
	    
	    target_v.z = molecules[i_m + Next_m].p.vz*dti;  // remove it after debugging
	    target_v.r = molecules[i_m + Next_m].p.vr*dti;  // remove it after debugging 
	    target_v.t = molecules[i_m + Next_m].p.vt*dti;  // remove it after debugging
	    
	    electrons[i_el + Next_el].p.vz -= M_m/(M_m + m_e)*v_rel.z;
	    electrons[i_el + Next_el].p.vr -= M_m/(M_m + m_e)*v_rel.r;
	    electrons[i_el + Next_el].p.vt -= M_m/(M_m + m_e)*v_rel.t;
	    
	    target_v.z += m_e/(M_m + m_e)*v_rel.z;
	    target_v.r += m_e/(M_m + m_e)*v_rel.r;
	    target_v.t += m_e/(M_m + m_e)*v_rel.t;
	    
	    /*                
	     *engcheck -= 0.5*m_e*(SQU(electrons[i_el + Next_el].p.vz) + SQU(electrons[i_el + Next_el].p.vr)+ SQU(electrons[i_el + Next_el].p.vt));
	     *engcheck -= 0.5*M_m*(SQU(target_v.z) + SQU(target_v.r)+ SQU(target_v.t));
	     
	     (*momcheck).z -= m_e*electrons[i_el + Next_el].p.vz + M_m*target_v.z;
	     (*momcheck).r -= m_e*electrons[i_el + Next_el].p.vr + M_m*target_v.r;
	     (*momcheck).t -= m_e*electrons[i_el + Next_el].p.vt + M_m*target_v.t;
	    */
	    
	    
	  }   //  if (W2_0 >= E_th && RAND < n_m*W_0*Si)
	  
	}   //for ( i_el = 0; i_el < n_el; i_el++)
      }   //if ( n_m && n_el )
      
      
      Next_el += n_el;
      Next_m += n_m;
      
    } //END for (int jz=0; jz<nz; jz++) 
  } //END for (int jr=0; jr<nr; jr++) 
}



/**********************************************************************

 Function:     coll_el_neutrals(...)
 Action:       neutrals-electrons collisions with production of ions (e + n -> 2e + i )
               1. We do inelastic collision with neutral in which we lose E_th
               2. We do ellastic collision with electron.

 Varied variables:
      electrons[], ...

 knm 31.07.01
     12.11.01  Adding Super-Super Neutrals ...

***********************************************************************/

void coll_el_neutrals_2D( Particle  neutrals[], size_t *nn, size_t ordcount_ntrls[], double M_n,
			  Particle  electrons[], size_t *ne, size_t ordcount_el[],
			  Particle  ions[], size_t *ni, int nr, int nz, int NZ, Reaction React,
			  Vec3d *momcheck, double *engcheck ) {
  
  const double m_e =1.;
  
  double W2, W, IW, W2_0, W_0, IW_0, v_proj_zr, ivp, cos_theta, sin_theta, cos_beta, sin_beta;
  //double W0;
  double cos_phi, sin_phi, psi;
  
  double vz, vr, vt; //double vx, vy, vz;
  Vec3d v_rel, delta_v_rel, target_v;
  
  static std::vector<size_t> gone_n;      //Not nice for parallelization
  static std::vector<size_t> sub_n2coll;  // (neither over spieces- or cell)
  static std::vector<size_t> sub_el2coll;
  
  double E_th  = React.Eth;
  double Emin  = React.Emin;
  //double Emax  = React.Emax;
  double Estep = React.Estep;
  int Epoints = React.N;
  int Eind;
  double S_i = 0.;
  
  size_t inzdlcl = 0;  
  size_t Next_el = 0;
  size_t Next_n  = 0;
  size_t shift_n = 0;
  
  double dti  = 1./dt_ion;      // Scaling factor 4 ion velocity

  for (int jr=0; jr<nr; jr++) {
    for (int jz=0; jz<nz; jz++) {
      
      size_t n_n  = ordcount_ntrls[jr*NZ+jz];
      size_t n_el = ordcount_el[jr*NZ+jz];
      
      //Play it safe -- resize everything to
      // the maximum num particles in the cell
      gone_n.resize(MAX(n_n,n_el)); 
      sub_n2coll.resize(MAX(n_n,n_el));
      sub_el2coll.resize(MAX(n_n,n_el));
      
      if ( n_n && n_el ) {
	size_t Nn_real  = 0;
	
	for (size_t i=0; i<n_n; i++ ) {
	  gone_n[i]= 0;
	  Nn_real += neutrals[i + Next_n].p.m;
	  sub_n2coll[i] = neutrals[i + Next_n].p.m;
	}
	//  fprintf(stderr,"!!! Node %d CELL %d Nn= %d", Node, j, Nn_real );
	
	
	if ( n_el<Nn_real ) {
	  size_t N2coll = Nn_real;
	  for (size_t i_el=0; i_el<n_el; i_el++ ) {
	    size_t k = (size_t)(N2coll*RAND);
	    if (k == N2coll) k--;   // index of "real particles"
	    
	    size_t i_n = 0;             
	    size_t i_next   = sub_n2coll[i_n];
	    
	    while ( i_next <= k ) {
	      i_next += sub_n2coll[++i_n];    
	    }
	    
	    --sub_n2coll[i_n];  
	    N2coll--;
	    
	    
	    // relative velocity
	    v_rel.z = electrons[i_el + Next_el].p.vz - neutrals[i_n + Next_n].p.vz*dti;   
	    v_rel.r = electrons[i_el + Next_el].p.vr - neutrals[i_n + Next_n].p.vr*dti;  
	    v_rel.t = electrons[i_el + Next_el].p.vt - neutrals[i_n + Next_n].p.vt*dti;
	    
	    // before scattering
	    W2_0 = SQU(v_rel.z) + SQU(v_rel.r) + SQU(v_rel.t);                
	    W_0  = sqrt(W2_0);
	    
	    // Linear fit for Cross-section
	    if ( W2_0 >= E_th ) {
	      Eind = (int)((W2_0 - Emin)/Estep);
	      if ( Eind < 0)
		S_i = React.CS[0];
	      else if ( Eind >= Epoints)
		S_i = React.CS[Epoints];
	      else
		S_i = React.CS[Eind] + (React.CS[Eind +1] - React.CS[Eind])*((W2_0 - Emin)/Estep - Eind);
	    }
	    S_i *= Ampl;
	    
	    // 1-exp(-n_n*W_0*Si*dt_coll)
	    if (W2_0 >= E_th && RAND < Nn_real*W_0*S_i/(2*jr+1) )  {// correct Si with 2D volume factor (2j+1) here
	      
	      *engcheck += 0.5*m_e*(SQU(electrons[i_el + Next_el].p.vz) + SQU(electrons[i_el + Next_el].p.vr) + SQU(electrons[i_el + Next_el].p.vt));
	      *engcheck += 0.5*M_n*(SQU(neutrals[i_n + Next_n].p.vz*dti) + SQU(neutrals[i_n + Next_n].p.vr*dti) + SQU(neutrals[i_n + Next_n].p.vt*dti)) - 0.5*m_e*M_n/(m_e+M_n)*E_th;
	      
	      (*momcheck).z += m_e*electrons[i_el + Next_el].p.vz + M_n*neutrals[i_n + Next_n].p.vz*dti;
	      (*momcheck).r += m_e*electrons[i_el + Next_el].p.vr + M_n*neutrals[i_n + Next_n].p.vr*dti;
	      (*momcheck).t += m_e*electrons[i_el + Next_el].p.vt + M_n*neutrals[i_n + Next_n].p.vt*dti;
	      
	      IW_0 = 1./W_0;
	      
	      W2 = W2_0 - E_th;
	      W  = sqrt(W2);                
	      
	      // After energy loss
	      delta_v_rel.z = v_rel.z*(W*IW_0 -1.);
	      delta_v_rel.r = v_rel.r*(W*IW_0 -1.);  
	      delta_v_rel.t = v_rel.t*(W*IW_0 -1.);
	      
	      v_rel.z += delta_v_rel.z;
	      v_rel.r += delta_v_rel.r;
	      v_rel.t += delta_v_rel.t;
	      
	      W  = MAX( W, 1.e-10);                         
	      IW = 1./W;
	      
	      v_proj_zr = sqrt( SQU(v_rel.z) + SQU(v_rel.r) );    
	      v_proj_zr = MAX( v_proj_zr, 1.e-10);              
	      ivp = 1./v_proj_zr;
	      
	      cos_theta = v_rel.t * IW;   
	      sin_theta = v_proj_zr*IW;
	      cos_beta = v_rel.z*ivp;     
	      sin_beta = v_rel.r*ivp;
	      
	      /* Sampling scattering angles in C.M. system */
	      // just 4 test!! Use = (2.+E-2.*pow((1+E), RAND))/E instead
	      cos_phi = 2.*RAND -1. ; 
	      sin_phi = sqrt(1. - SQU(cos_phi));
	      psi  = TWOPI*RAND;
	      
	      vt = W*cos_phi;              
	      vz = W*sin_phi*cos(psi);
	      vr = W*sin_phi*sin(psi);
	      
	      v_rel.z -= vz*cos_beta*cos_theta - vr*sin_beta + vt*cos_beta*sin_theta;  
	      v_rel.r -= vz*sin_beta*cos_theta + vr*cos_beta + vt*sin_beta*sin_theta;  
	      v_rel.t -= vt*cos_theta - vz*sin_theta;
	      
	      v_rel.z -= delta_v_rel.z;
	      v_rel.r -= delta_v_rel.r;
	      v_rel.t -= delta_v_rel.t;
	      
	      target_v.z = neutrals[i_n + Next_n].p.vz*dti;
	      target_v.r = neutrals[i_n + Next_n].p.vr*dti;
	      target_v.t = neutrals[i_n + Next_n].p.vt*dti;
	      
	      electrons[i_el + Next_el].p.vz -= M_n/(M_n + m_e)*v_rel.z;
	      electrons[i_el + Next_el].p.vr -= M_n/(M_n + m_e)*v_rel.r;
	      electrons[i_el + Next_el].p.vt -= M_n/(M_n + m_e)*v_rel.t;
	      
	      target_v.z += m_e/(M_n + m_e)*v_rel.z; 
	      target_v.r += m_e/(M_n + m_e)*v_rel.r;    
	      target_v.t += m_e/(M_n + m_e)*v_rel.t;    

	      if ( ((*ni) + inzdlcl) >= NPART) {
		printf("Error in coll_el_neutrals_2D: Particle array overflow (ions, case1)\n");
		exit(1);
	      }

	      // 2D: add also z-component of particle position
	      ions[*ni + inzdlcl].p.z = neutrals[i_n + Next_n].p.z;
	      ions[*ni + inzdlcl].p.r = neutrals[i_n + Next_n].p.r;  
	      
	      // scale it back
	      ions[*ni + inzdlcl].p.vz = target_v.z*dt_ion;
	      ions[*ni + inzdlcl].p.vr = target_v.r*dt_ion;  
	      ions[*ni + inzdlcl].p.vt = target_v.t*dt_ion;
	      ions[*ni + inzdlcl].p.m = 1;             
	      
	      *engcheck -= 0.5*(M_n-m_e)*(SQU(ions[*ni + inzdlcl].p.vz*dti) + SQU(ions[*ni + inzdlcl].p.vr*dti) + SQU(ions[*ni + inzdlcl].p.vt*dti));
	      
	      (*momcheck).z -= (M_n-m_e)*ions[*ni + inzdlcl].p.vz*dti;
	      (*momcheck).r -= (M_n-m_e)*ions[*ni + inzdlcl].p.vr*dti;
	      (*momcheck).t -= (M_n-m_e)*ions[*ni + inzdlcl].p.vt*dti;
	      
	      
	      /* electron-neutral inelastic collision with a loss of E_th energy
	       * and birth of the ion  is done !
	       * now we'll do electron-electron elastic collision
	       */
	      
	      v_rel.z = electrons[i_el + Next_el].p.vz -  target_v.z;
	      v_rel.r = electrons[i_el + Next_el].p.vr -  target_v.r;    
	      v_rel.t = electrons[i_el + Next_el].p.vt -  target_v.t;
	      
	      
	      W2 = SQU(v_rel.z) + SQU(v_rel.r) + SQU(v_rel.t);
	      W  = sqrt(W2);                     
	      W  = MAX( W, 1.e-10);            
	      IW = 1./W;
	      
	      v_proj_zr = sqrt( SQU(v_rel.z) + SQU(v_rel.r) );  
	      v_proj_zr = MAX( v_proj_zr, 0.000001);              
	      ivp       = 1./v_proj_zr;
	      
	      cos_theta = v_rel.t * IW;         
	      sin_theta = v_proj_zr*IW;
	      cos_beta = v_rel.z*ivp;     
	      sin_beta = v_rel.r*ivp;
	      
	      cos_phi = 2.*RAND -1. ;
	      sin_phi = sqrt(1 - SQU(cos_phi));
	      psi  = TWOPI*RAND;
	      
	      vt = W*cos_phi;           
	      vz = W*sin_phi*cos(psi);
	      vr = W*sin_phi*sin(psi);
	      
	      v_rel.z -= vz*cos_beta*cos_theta - vr*sin_beta + vt*cos_beta*sin_theta;   
	      v_rel.r -= vz*sin_beta*cos_theta + vr*cos_beta + vt*sin_beta*sin_theta;   
	      v_rel.t -= vt*cos_theta - vz*sin_theta;
	      
	      // Note inversed sign due inversed delta be4
	      electrons[i_el + Next_el].p.vz -= 0.5*v_rel.z; 
	      electrons[i_el + Next_el].p.vr -= 0.5*v_rel.r;          
	      electrons[i_el + Next_el].p.vt -= 0.5*v_rel.t;
	      
	      target_v.z += 0.5*v_rel.z;
	      target_v.r += 0.5*v_rel.r;
	      target_v.t += 0.5*v_rel.t;

	      if ( ((*ne) + inzdlcl) >= NPART) {
		printf("Error in coll_el_neutrals_2D: Particle array overflow (electrons, case1)\n");
		exit(1);
	      }	      
	      
	      electrons[*ne + inzdlcl].p.z = neutrals[i_n + Next_n].p.z;
	      electrons[*ne + inzdlcl].p.r = neutrals[i_n + Next_n].p.r;
	      electrons[*ne + inzdlcl].p.vz = target_v.z;
	      electrons[*ne + inzdlcl].p.vr = target_v.r;
	      electrons[*ne + inzdlcl].p.vt = target_v.t;
	      electrons[*ne + inzdlcl].p.m = 1;      
	      
	      
	      *engcheck -= 0.5*m_e*(SQU(electrons[i_el + Next_el].p.vz) + SQU(electrons[i_el + Next_el].p.vr) + SQU(electrons[i_el + Next_el].p.vt));
	      *engcheck -= 0.5*m_e*(SQU(electrons[*ne + inzdlcl].p.vz) + SQU(electrons[*ne + inzdlcl].p.vr) + SQU(electrons[*ne + inzdlcl].p.vt));
	      
	      (*momcheck).z -= m_e*(electrons[i_el + Next_el].p.vz + electrons[*ne + inzdlcl].p.vz);
	      (*momcheck).r -= m_e*(electrons[i_el + Next_el].p.vr + electrons[*ne + inzdlcl].p.vr);
	      (*momcheck).t -= m_e*(electrons[i_el + Next_el].p.vt + electrons[*ne + inzdlcl].p.vt);
	      
	      
	      inzdlcl++;
	      
	      
	      if (! --(neutrals[i_n + Next_n].p.m) ) {
		gone_n[i_n]= 1;      
	      }
	      
	    }   //  if (W2_0 >= E_th && RAND < n_n*W_0*Si)
	  }  //for ( i_el = 0; i_el < n_el; i_el++)
          
          
	  for (size_t i = 0; i < n_n; i++) {
	    if ( gone_n[i])  shift_n++ ;
	    else neutrals[Next_n + i - shift_n] = neutrals[Next_n + i]; 
	  }
	  
	}   // if (n_el < n_n)
	
	/* Case when n_n < n_el*/
	else {
	  size_t N2coll = n_el;
	  for (size_t i=0; i<n_el; i++ ) {	                
	    sub_el2coll[i] = 1; 
	  }
          
	  for (size_t i_n=0; i_n<n_n; i_n++)
	    for (size_t i_n_sub=0; i_n_sub<sub_n2coll[i_n]; i_n_sub++) {
	      
	      size_t k = (size_t)(N2coll*RAND);
	      if ( k==N2coll ) k--;          
	      
	      size_t i_el = 0;
	      size_t i_next = sub_el2coll[i_el];
	      
	      while ( i_next<=k) {
		i_next += sub_el2coll[++i_el];
	      }
	      
	      --sub_el2coll[i_el];   
	      N2coll--;
	      
	      v_rel.z= electrons[i_el + Next_el].p.vz - neutrals[i_n + Next_n].p.vz*dti;
	      v_rel.r= electrons[i_el + Next_el].p.vr - neutrals[i_n + Next_n].p.vr*dti;   
	      v_rel.t= electrons[i_el + Next_el].p.vt - neutrals[i_n + Next_n].p.vt*dti;
	      
	      W2_0 = SQU(v_rel.z) + SQU(v_rel.r) + SQU(v_rel.t);  
	      W_0  = sqrt(W2_0);
	      
	      
	      // Linear fit for Cross-section
	      if ( W2_0 >= E_th ) {
		Eind = (int)((W2_0 - Emin)/Estep);
		if ( Eind < 0)
		  S_i = React.CS[0];
		else if ( Eind >= Epoints)
		  S_i = React.CS[Epoints];
		else
		  S_i = React.CS[Eind] + (React.CS[Eind +1] - React.CS[Eind])*((W2_0 - Emin)/Estep - Eind);
	      }
	      S_i *= Ampl;
	      
	      if (W2_0 >= E_th && RAND < n_el*W_0*S_i/(2*jr+1) ) { // correct Si with 2D volume factor (2j+1) here
		
		
		*engcheck += 0.5*m_e*(SQU(electrons[i_el + Next_el].p.vz) + SQU(electrons[i_el + Next_el].p.vr) + SQU(electrons[i_el + Next_el].p.vt));
		*engcheck += 0.5*M_n*(SQU(neutrals[i_n + Next_n].p.vz*dti) + SQU(neutrals[i_n + Next_n].p.vr*dti) + SQU(neutrals[i_n + Next_n].p.vt*dti)) - 0.5*m_e*M_n/(m_e+M_n)*E_th;
		
		(*momcheck).z += m_e*electrons[i_el + Next_el].p.vz + M_n*neutrals[i_n + Next_n].p.vz*dti;
		(*momcheck).r += m_e*electrons[i_el + Next_el].p.vr + M_n*neutrals[i_n + Next_n].p.vr*dti;
		(*momcheck).t += m_e*electrons[i_el + Next_el].p.vt + M_n*neutrals[i_n + Next_n].p.vt*dti;
		
		IW_0 = 1./W_0;
		
		W2 = W2_0 - E_th;
		W  = sqrt(W2);                
		  
		// After energy loss
		delta_v_rel.z = v_rel.z*(W*IW_0 -1.);
		delta_v_rel.r = v_rel.r*(W*IW_0 -1.);       
		delta_v_rel.t = v_rel.t*(W*IW_0 -1.);
		
		v_rel.z += delta_v_rel.z;
		v_rel.r += delta_v_rel.r;
		v_rel.t += delta_v_rel.t;
		
		W  = MAX( W, 1.e-10);   
		IW = 1./W;
		
		
		v_proj_zr = sqrt( SQU(v_rel.z) + SQU(v_rel.r) );    
		v_proj_zr = MAX( v_proj_zr, 1.e-10);              
		ivp       = 1./v_proj_zr;
		
		cos_theta = v_rel.t * IW;             
		sin_theta = v_proj_zr*IW;
		cos_beta = v_rel.z*ivp;         
		sin_beta = v_rel.r*ivp;
		
		/* Sampling scattering angles in C.M. system */
		// just 4 test!! Use = (2.+E-2.*pow((1+E), RAND))/E instead
		cos_phi = 2.*RAND -1. ; 
		sin_phi = sqrt(1. - SQU(cos_phi));
		psi  = TWOPI*RAND;
		
		vt = W*cos_phi;              
		vz = W*sin_phi*cos(psi);
		vr = W*sin_phi*sin(psi);
		
		v_rel.z -= vz*cos_beta*cos_theta - vr*sin_beta + vt*cos_beta*sin_theta; 
		v_rel.r -= vz*sin_beta*cos_theta + vr*cos_beta + vt*sin_beta*sin_theta;  
		v_rel.t -= vt*cos_theta - vz*sin_theta;
		
		v_rel.z -= delta_v_rel.z;
		v_rel.r -= delta_v_rel.r;
		v_rel.t -= delta_v_rel.t;
		
		target_v.z = neutrals[i_n + Next_n].p.vz*dti;
		target_v.r = neutrals[i_n + Next_n].p.vr*dti;
		target_v.t = neutrals[i_n + Next_n].p.vt*dti;
		
		electrons[i_el + Next_el].p.vz -= M_n/(M_n + m_e)*v_rel.z;
		electrons[i_el + Next_el].p.vr -= M_n/(M_n + m_e)*v_rel.r;
		electrons[i_el + Next_el].p.vt -= M_n/(M_n + m_e)*v_rel.t;
		
		target_v.z += m_e/(M_n + m_e)*v_rel.z;  
		target_v.r += m_e/(M_n + m_e)*v_rel.r;    
		target_v.t += m_e/(M_n + m_e)*v_rel.t;    
		
		if ( ((*ni) + inzdlcl) >= NPART) {
		  printf("Error in coll_el_neutrals_2D: Particle array overflow (ions, case2)\n");
		  exit(1);
		}
		
		ions[*ni + inzdlcl].p.z = neutrals[i_n + Next_n].p.z;
		ions[*ni + inzdlcl].p.r = neutrals[i_n + Next_n].p.r; 
		ions[*ni + inzdlcl].p.vz = target_v.z*dt_ion;
		ions[*ni + inzdlcl].p.vr = target_v.r*dt_ion; 
		ions[*ni + inzdlcl].p.vt = target_v.t*dt_ion;
		ions[*ni + inzdlcl].p.m = 1;                             
		
		*engcheck -= 0.5*(M_n-m_e)*(SQU(ions[*ni + inzdlcl].p.vz*dti) + SQU(ions[*ni + inzdlcl].p.vr*dti) + SQU(ions[*ni + inzdlcl].p.vt*dti));
		
		(*momcheck).z -= (M_n-m_e)*ions[*ni + inzdlcl].p.vz*dti;
		(*momcheck).r -= (M_n-m_e)*ions[*ni + inzdlcl].p.vr*dti;
		(*momcheck).t -= (M_n-m_e)*ions[*ni + inzdlcl].p.vt*dti;
		
		/*
		 * electron-neutral inelastic collision with a loss of E_th energy
		 * and birth of the ion  is done !
		 * now we'll do electron-electron elastic collission
		 */
		
		v_rel.z= electrons[i_el + Next_el].p.vz -  target_v.z;
		v_rel.r= electrons[i_el + Next_el].p.vr -  target_v.r; 
		v_rel.t= electrons[i_el + Next_el].p.vt -  target_v.t;
		
		W2 = SQU(v_rel.z) + SQU(v_rel.r) + SQU(v_rel.t);
		W  = sqrt(W2);                      
		W  = MAX( W, 1.e-10);             
		IW = 1./W;
		
		v_proj_zr = sqrt( SQU(v_rel.z) + SQU(v_rel.r) );  
		v_proj_zr = MAX( v_proj_zr, 0.000001);              
		ivp       = 1./v_proj_zr;
		
		cos_theta = v_rel.t * IW;    
		sin_theta = v_proj_zr*IW;
		cos_beta = v_rel.z*ivp;       
		sin_beta = v_rel.r*ivp;
		
		cos_phi = 2.*RAND -1. ;
		sin_phi = sqrt(1 - SQU(cos_phi));
		psi  = TWOPI*RAND;
		
		vt = W*cos_phi;             
		vz = W*sin_phi*cos(psi);
		vr = W*sin_phi*sin(psi);
		
		v_rel.z -= vz*cos_beta*cos_theta - vr*sin_beta + vt*cos_beta*sin_theta; 
		v_rel.r -= vz*sin_beta*cos_theta + vr*cos_beta + vt*sin_beta*sin_theta;  
		v_rel.t -= vt*cos_theta - vz*sin_theta;
		
		electrons[i_el + Next_el].p.vz -= 0.5*v_rel.z;
		electrons[i_el + Next_el].p.vr -= 0.5*v_rel.r;  
		electrons[i_el + Next_el].p.vt -= 0.5*v_rel.t;
		
		target_v.z += 0.5*v_rel.z;
		target_v.r += 0.5*v_rel.r;
		target_v.t += 0.5*v_rel.t;

		if ( ((*ne) + inzdlcl) >= NPART) {
		  printf("Error in coll_el_neutrals_2D: Particle array overflow (electrons, case2)\n");
		  exit(1);
		}		
		
		electrons[*ne + inzdlcl].p.z = neutrals[i_n + Next_n].p.z;
		electrons[*ne + inzdlcl].p.r = neutrals[i_n + Next_n].p.r;
		electrons[*ne + inzdlcl].p.vz = target_v.z;
		electrons[*ne + inzdlcl].p.vr = target_v.r;
		electrons[*ne + inzdlcl].p.vt = target_v.t;
		electrons[*ne + inzdlcl].p.m  = 1;
		
		
		*engcheck -= 0.5*m_e*(SQU(electrons[i_el + Next_el].p.vz) + SQU(electrons[i_el + Next_el].p.vr) + SQU(electrons[i_el + Next_el].p.vt));
		*engcheck -= 0.5*m_e*(SQU(electrons[*ne + inzdlcl].p.vz) + SQU(electrons[*ne + inzdlcl].p.vr) + SQU(electrons[*ne + inzdlcl].p.vt));
		
		(*momcheck).z -= m_e*(electrons[i_el + Next_el].p.vz + electrons[*ne + inzdlcl].p.vz);
		(*momcheck).r -= m_e*(electrons[i_el + Next_el].p.vr + electrons[*ne + inzdlcl].p.vr);
		(*momcheck).t -= m_e*(electrons[i_el + Next_el].p.vt + electrons[*ne + inzdlcl].p.vt);
		
		
		inzdlcl++;
		if (! --neutrals[i_n + Next_n].p.m ) {
		  gone_n[i_n]= 1; 
		}
		
	      }   //  if (W2_0 >= E_th && RAND < n_n*W_0*Si)
	    }  //for ( i_n_sub = 0; i_n_sub < sub_n2coll[i_n]; i_n_sub++)
	  
	  
	  for (size_t i=0; i<n_n; i++ ) {
	    if ( gone_n[i] )  shift_n++ ;
	    else neutrals[Next_n + i - shift_n] = neutrals[Next_n + i]; 
	  }
	  
	}  // else if n_el < n_n
      }  //if ( n_n && n_el )
      else {
	for (size_t i=0; i<n_n; i++ ) {
	  neutrals[Next_n + i - shift_n] = neutrals[Next_n + i];
	}
      }
      
      Next_el += n_el;
      Next_n += n_n;
    }    
  }  //for (j=0; j < Ngrid;  j++)
  
  *nn -= shift_n;
  *ne += inzdlcl;  
  *ni += inzdlcl;  
  
}
