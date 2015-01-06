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

  push.cpp:
  Particle mover; eletrostatic treatment only (optional B_ext)

***********************************************************************/

#include <stdio.h>
#include <math.h>

#include  "pic.h"
#include  "dim.h"
#include  "mydef.h" 

void  push_2D( Particle pa[], const double Eg_r[], const double Eg_z[], size_t np, int NZ ) {
  int       j, k;
  double    hr,  Er;
  double    hz,  Ez;
  double    r0, ca, sa, vr;
  
  for(size_t n=0; n < np; n++) {
    hr  = pa[n].p.r;
    j   = (int)hr;
    hr -= j;       
    
    hz  = pa[n].p.z;
    k   = (int)hz;
    hz -= k;   
    
    // Interpolation
    Ez  = (1-hr)*(1-hz)*Eg_z[j*NZ+k] + (1-hr)*hz*Eg_z[j*NZ+(k+1)] + hr*(1-hz)*Eg_z[(j+1)*NZ+k] + hr*hz*Eg_z[(j+1)*NZ+(k+1)]; 
    Er  = (1-hr)*(1-hz)*Eg_r[j*NZ+k] + (1-hr)*hz*Eg_r[j*NZ+(k+1)] + hr*(1-hz)*Eg_r[(j+1)*NZ+k] + hr*hz*Eg_r[(j+1)*NZ+(k+1)];
    
    // Acceleration
    pa[n].p.vz -= 2.*Ez;
    pa[n].p.vr -= 2.*Er;
    
    // Move particle 
    r0 = pa[n].p.r;
    pa[n].p.z += pa[n].p.vz; 
    pa[n].p.r = sqrt( SQU(r0 + pa[n].p.vr) + SQU(pa[n].p.vt) );
    
    // Rotate coordinate system
    if ( pa[n].p.r > 1.e-8 ) {
      sa = pa[n].p.vt/pa[n].p.r;
      ca = (r0 + pa[n].p.vr)/pa[n].p.r;
    }
    else {
      if ( pa[n].p.vr < -1.e-8 ) {
	sa = 0.;
	ca = -1.;
      }
      else {
	sa = 0.;
	ca = 1.;
      }
    }
    
    vr = pa[n].p.vr;
    pa[n].p.vr = ca * vr + sa * pa[n].p.vt;
    pa[n].p.vt = -sa * vr + ca * pa[n].p.vt;    
  }
}

void  push_magnetic_onlyBz_2D( Particle pa[], const double Eg_r[], const double Eg_z[], 
			       const double Bz, double factor, size_t np, int NZ ) {
  
  int       j, k;
  double    hr,  Er;
  double    hz,  Ez;
  double    r0, ca, sa, vr;
  
  double B = factor*Bz;
  double B_alt = 2*B/(1+B*B); 
  double vza, vra;
  double vrb, vtb;
  double vrc;
  
  for(size_t n=0; n < np; n++) {
    hr  = pa[n].p.r;
    j   = (int)hr;
    hr -= j;       
    
    hz  = pa[n].p.z;
    k   = (int)hz;
    hz -= k;   
    
    // Interpolation
    Ez  = (1-hr)*(1-hz)*Eg_z[j*NZ+k] + (1-hr)*hz*Eg_z[j*NZ+(k+1)] + hr*(1-hz)*Eg_z[(j+1)*NZ+k] + hr*hz*Eg_z[(j+1)*NZ+(k+1)]; 
    Er  = (1-hr)*(1-hz)*Eg_r[j*NZ+k] + (1-hr)*hz*Eg_r[j*NZ+(k+1)] + hr*(1-hz)*Eg_r[(j+1)*NZ+k] + hr*hz*Eg_r[(j+1)*NZ+(k+1)];
    
    // Acceleration: Boris
    vza = pa[n].p.vz - Ez;	  
    vra = pa[n].p.vr - Er;	  
    
    vrb = vra - pa[n].p.vt*B;
    vtb = pa[n].p.vt + vra*B;
    
    vrc = vra - vtb*B_alt;
    
    pa[n].p.vt += vrb*B_alt;
    pa[n].p.vz = vza - Ez;
    pa[n].p.vr = vrc - Er;
    
    // Move particle 
    r0 = pa[n].p.r;
    pa[n].p.z += pa[n].p.vz; 
    pa[n].p.r = sqrt( SQU(r0 + pa[n].p.vr) + SQU(pa[n].p.vt) );
    
    // Rotate coordinate system
    if ( pa[n].p.r > 1.e-8 ) {
      sa = pa[n].p.vt/pa[n].p.r;
      ca = (r0 + pa[n].p.vr)/pa[n].p.r;
    }
    else {
      if ( pa[n].p.vr < -1.e-8 ) {
	sa = 0.;
	ca = -1.;
      }
      else {
	sa = 0.;
	ca = 1.;
      }
    }   
    vr = pa[n].p.vr;
    pa[n].p.vr = ca * vr + sa * pa[n].p.vt;
    pa[n].p.vt = -sa * vr + ca * pa[n].p.vt;   
  }

}

void  push_magnetic_2D( Particle pa[], const double Eg_r[], const double Eg_z[], 
			const double Bextz, const double Bextt, double factor, size_t np, int NZ ) {

  int       j, k;
  double    hr,  Er;
  double    hz,  Ez;
  double    r0, ca, sa, vr;
  
  double Bz = factor*Bextz;
  double Bt = factor*Bextt;
  double B_alt = 2./(1.+Bz*Bz+Bt*Bt); 
  double vza, vra;
  double vzb, vrb, vtb;
  double vzc, vrc;
  
  for(size_t n=0; n < np; n++) {
    hr  = pa[n].p.r;
    j   = (int)hr;
    hr -= j;       
    
    hz  = pa[n].p.z;
    k   = (int)hz;
    hz -= k;   
    
    // Interpolation
    Ez  = (1-hr)*(1-hz)*Eg_z[j*NZ+k] + (1-hr)*hz*Eg_z[j*NZ+(k+1)] + hr*(1-hz)*Eg_z[(j+1)*NZ+k] + hr*hz*Eg_z[(j+1)*NZ+(k+1)]; 
    Er  = (1-hr)*(1-hz)*Eg_r[j*NZ+k] + (1-hr)*hz*Eg_r[j*NZ+(k+1)] + hr*(1-hz)*Eg_r[(j+1)*NZ+k] + hr*hz*Eg_r[(j+1)*NZ+(k+1)];
    
    // Acceleration: Boris
    vza = pa[n].p.vz - Ez;	  
    vra = pa[n].p.vr - Er;	  
    
    vzb = vza - vra*Bt;
    vrb = vra - pa[n].p.vt*Bz + vza*Bt;
    vtb = pa[n].p.vt + vra*Bz;
    
    vzc = vza - vrb*Bt*B_alt;
    vrc = vra - (vtb*Bz - vzb*Bt)*B_alt;
    pa[n].p.vt += vrb*Bz*B_alt;
    
    pa[n].p.vz = vzc - Ez;
    pa[n].p.vr = vrc - Er;
    
    // Move particle 
    r0 = pa[n].p.r;
    pa[n].p.z += pa[n].p.vz; 
    pa[n].p.r = sqrt( SQU(r0 + pa[n].p.vr) + SQU(pa[n].p.vt) );  
    
    // Rotate coordinate system
    if ( pa[n].p.r > 1.e-8 ) {
      sa = pa[n].p.vt/pa[n].p.r;
      ca = (r0 + pa[n].p.vr)/pa[n].p.r;
    }
    else {
      if ( pa[n].p.vr < -1.e-8 ) {
	sa = 0.;
	ca = -1.;
      }
      else {
	sa = 0.;
	ca = 1.;
      }
    } 
    vr = pa[n].p.vr;
    pa[n].p.vr = ca * vr + sa * pa[n].p.vt;
    pa[n].p.vt = -sa * vr + ca * pa[n].p.vt;
  }
}

void  push_neutrals_2D( Particle pa[], size_t np ) {
  
  double    r0, ca, sa, vr;
  
  for(size_t n=0; n < np; n++) {
    // No acceleration
    // Move particle 
    r0 = pa[n].p.r;
    pa[n].p.z += pa[n].p.vz; 
    pa[n].p.r = sqrt( SQU(r0 + pa[n].p.vr) + SQU(pa[n].p.vt) );
    
    // Rotate coordinate system
    if ( pa[n].p.r > 1.e-8 ) {
      sa = pa[n].p.vt/pa[n].p.r;
      ca = (r0 + pa[n].p.vr)/pa[n].p.r;
    }
    else {
      if ( pa[n].p.vr < -1.e-8 ) {
	sa = 0.;
	ca = -1.;
      }
      else {
	sa = 0.;
	ca = 1.;
      }
    }
    vr = pa[n].p.vr;
    pa[n].p.vr = ca * vr + sa * pa[n].p.vt;
    pa[n].p.vt = -sa * vr + ca * pa[n].p.vt;   
  }

}
