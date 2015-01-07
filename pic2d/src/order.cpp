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

  order.cpp:
  Orders particles according to cells, after Konstantin Matyash

***********************************************************************/

#include  "pic.h"
#include  "dim.h"

#include <stdio.h>
#include <stdlib.h>

#include <vector>

void order_2D(Particle pa[], size_t np, size_t ordcount[], int nr, int nz, int NZ) {
  //Static storage list
  // Note: Static not nice for parallelization...
  static std::vector<Particle>* temP = NULL;
  static int nr_def(-1), nz_def(-1);
  if (!temP) {
    //Initialize at first pass
    nr_def=nr; nz_def=nz;
    temP = new std::vector<Particle>[nr*nz];
  }
  //Sanity check
  if( nr!=nr_def or nz!=nz_def ) {
    fprintf(stderr,"nr or nz changed!");
    fflush(stderr);
    printf("Error in order_2D(): nr or nz changed between calls, temP invalid!\n");
    exit(1);
  }
  
  // Initialize ordcount
  for (int j=0; j<nr; j++) {
    for (int k=0; k<nz; k++){
      ordcount[j*NZ+k]= 0;
    }
  }
  
  //Sort particles
  int ir, iz;
  for (size_t n=0; n<np; n++) {
    ir = (int)pa[n].p.r;
    iz = (int)pa[n].p.z;

    //Sanity check -- we should be inside the grid   
    if ( (ir<0) || (ir>=nr) || (iz<0) || (iz>=nz) ) {
      //Edge cases
      if (pa[n].p.r == (double) nr) ir--;
      else if (pa[n].p.z == (double) nz) iz--;
      //Just wrong!
      else {
	fprintf(stderr," !!!ERROR!!! %d %d\n", ir, iz);
	fflush(stderr);

	printf("Error in order_2D(): particle position out-of-range, r=%g, z=%g, (ir,iz)= (%d,%d)\n", pa[n].p.r, pa[n].p.z, ir, iz);
	
	exit(1);
      }
    }

    temP[ir*nz+iz].push_back(pa[n]);
    ordcount[ir*NZ+iz]++;
  }
  
  //Stuff them back into pa, but now in the correct order
  size_t i=0;
  for (int cidx = 0; cidx < nr*nz; cidx++) {
    for (size_t l = 0; l < temP[cidx].size(); l++) {
      pa[i++] = temP[cidx][l];
    }
    //Clear the cell vector.
    // This will on *most* STL implementations
    // not lead to dealloc of memory!
    // http://www.informit.com/guides/content.aspx?g=cplusplus&seqNum=434 (for example)
    temP[cidx].clear();
  } 
  
}
