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

  moms.cpp:
  Calculates average macro-quantities (different moments)

***********************************************************************/

#include  "pic.h"
#include  "mydef.h"

#include <cstddef>

void aver_moments_2D( Moments mom[], int n_av,
		      Particle pa[], size_t np,
		      int nr, int nz, int NZ ) {
  
  if (n_av == 0) {
    for (int j=0; j<nr; j++) {
      for (int k=0; k<nz; k++) {
	mom[j*NZ+k].uz = mom[j*NZ+k].ur = mom[j*NZ+k].ut = 0.;
	mom[j*NZ+k].tz = mom[j*NZ+k].tr = mom[j*NZ+k].tt = 0.;
	mom[j*NZ+k].n = 0;
      }
    }
  }
  
  for (size_t n=0; n<np; n++) {
    int j = (int)pa[n].p.r;
    int k = (int)pa[n].p.z;
    
    double vz = pa[n].p.vz;
    double vr = pa[n].p.vr;
    double vt = pa[n].p.vt;
    mom[j*NZ+k].uz += vz;
    mom[j*NZ+k].ur += vr;
    mom[j*NZ+k].ut += vt;
    
    double vz2 = SQU(vz);
    double vr2 = SQU(vr);
    double vt2 = SQU(vt);
    mom[j*NZ+k].tz += vz2;
    mom[j*NZ+k].tr += vr2;
    mom[j*NZ+k].tt += vt2;
    
    mom[j*NZ+k].n++;
  }  
}

void aver_diagn_2D( const double dens[], double dens_av[],
		    Particle pa[], double temp_av[], double np_av[],
		    size_t np, int n_av, int nr, int nz, int NR, int NZ ) {
  
  if (n_av == 0) {
    for (int j=0; j<nr; j++) {
      for (int k=0; k<nz; k++) {
	temp_av[j*NZ+k] = 0.;
	np_av[j*NZ+k] = 0.;
      }
    }
    for (int j=0; j<NR; j++) {
      for (int k=0; k<NZ; k++) {
	dens_av[j*NZ+k] = dens[j*NZ+k];
      }
    }
  }
  
  else {
    for (int j=0; j<NR; j++) {
      for (int k=0; k<NZ; k++) {
	dens_av[j*NZ+k] += dens[j*NZ+k];
      }
    }
  }


  for (size_t n=0; n<np; n++) {
    int j = (int)pa[n].p.r;
    int k = (int)pa[n].p.z;

    double vz = pa[n].p.vz;
    double vr = pa[n].p.vr;
    double vt = pa[n].p.vt;
    
    temp_av[j*NZ+k] += SQU(vz) + SQU(vr) + SQU(vt);
    np_av[j*NZ+k]++;
  }
}




void aver_moments_SN_2D( Moments mom[], int n_av,
			 Particle pa[], size_t np,
			 int nr, int nz, int NZ ) {  
  
  if (n_av == 0) {
    for (int j=0; j<nr; j++) {
      for (int k=0; k<nz; k++) {
	mom[j*NZ+k].uz = mom[j*NZ+k].ur = mom[j*NZ+k].ut = 0.;
	mom[j*NZ+k].tz = mom[j*NZ+k].tr = mom[j*NZ+k].tt = 0.;
	mom[j*NZ+k].n = 0;
      }
    }
  }

  for (size_t n=0; n<np; n++) {
    int j = (int)pa[n].p.r;
    int k = (int)pa[n].p.z;
    int real_ones  = pa[n].p.m;
    
    double vz = pa[n].p.vz;
    double vr = pa[n].p.vr;
    double vt = pa[n].p.vt;
    mom[j*NZ+k].uz += vz*real_ones;
    mom[j*NZ+k].ur += vr*real_ones;
    mom[j*NZ+k].ut += vt*real_ones;
    
    double vz2 = SQU(vz);
    double vr2 = SQU(vr);
    double vt2 = SQU(vt);
    mom[j*NZ+k].tz += vz2*real_ones;
    mom[j*NZ+k].tr += vr2*real_ones;
    mom[j*NZ+k].tt += vt2*real_ones;
    
    mom[j*NZ+k].n += real_ones;
  }
  //don't forget to scale mean values on Super particles
}
