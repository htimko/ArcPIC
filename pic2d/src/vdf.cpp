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

  vdf.cpp:
  Calculates the velocity distribution function

***********************************************************************/

#include <cmath>
#include <cstddef>

#include  "pic.h"
#include  "mydef.h"
#include  "dim.h"



void vel_dst_along_2D( double* fvz, double* fvr, double* fvabs, Particle pa[], 
		       size_t np, int nr, int nz, double dvt, int n_av ) {
  
  int nvall = nr*nz*Nvdst/4;  
  int offset = Nvdst/2; // symmetric around v=0
  
  if (n_av == 0) {
    for (int m=0; m<nvall; m++) {
      fvz[m] = fvr[m] = fvabs[m] = 0;
    }
  }

  for (size_t n=0; n<np; n++) {
    int rr = (int)pa[n].p.r;
    int zz = (int)pa[n].p.z;
    
    int jj, kk;
    if ((jj = rr/2) == (rr+1)/2) { // every second grid
      if ((kk = zz/2) == (zz+1)/2) {
	
	int m = (int)(10*pa[n].p.vz*dvt+offset);
	if (m>=0 && m<Nvdst) fvz[m + Nvdst*(jj*(nz/2)+kk)] += 1.;

	m = (int)(10*pa[n].p.vr*dvt+offset);
	if (m>=0 && m<Nvdst) fvr[m + Nvdst*(jj*(nz/2)+kk)] += 1.;
	
	// absolute value distribution
	double v2 = SQU(pa[n].p.vz) + SQU(pa[n].p.vr) + SQU(pa[n].p.vt);
	double v = sqrt(v2);
	m = (int)(10*v*dvt+offset);
	if (m>=0 && m<Nvdst) fvabs[m + Nvdst*(jj*(nz/2)+kk)] += 1./MAX(v2, 1.e-10);
      }
    }
    
  }
  
}
