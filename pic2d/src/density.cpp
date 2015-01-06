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

  density.cpp:
  Density calculation, volume dependence in V_cell

***********************************************************************/

#include  "pic.h"

#include <cstddef>

void  density_2D( double dens[], int nr, int nz, int NR, int NZ, double V_cell[],
		  Particle pa[], double qp, size_t np, int sort, int Sorts ) {
  
  // Initialise
  for (int j=0; j<NR; j++) {
    for (int k=0; k<NZ; k++) {	  
      dens[j*NZ+k] = 0.; 
    }
  }


  // Positive ions: only this
  for (size_t n=0; n<np; n++) {
    double hr = pa[n].p.r;
    int j  = (int)hr;
    hr -= j;
    
    double hz = pa[n].p.z;
    int k  = (int)hz;
    hz -= k;
    
    dens[j*NZ+k]         += (1-hr)*(1-hz)/V_cell[j];
    dens[(j+1)*NZ+k]     += hr*(1-hz)/V_cell[j+1];
    dens[j*NZ+(k+1)]     += (1-hr)*hz/V_cell[j];
    dens[(j+1)*NZ+(k+1)] += hr*hz/V_cell[j+1];
  }
  
  /*  Multipying by particles charge	   */
  if (sort == Sorts-1 ) {
    for(int j=0; j<NR*NZ; j++)  dens[j] *= qp;
  }
}
