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

  efield.cpp:
  Electric field calcuation

***********************************************************************/

#include  "dim.h"

#define XTRN extern
#include  "var.h"
#undef XTRN


void electric_field_2D( double Phi[], double E_grid_r[], double E_grid_z[], 
			double E_ion_r[], double E_ion_z[], 
			int nr, int nz, int NR, int NZ ) {

  for (int j=0; j<NR; j++) {
    for (int k=0; k<NZ; k++) {

      if ( BC == 0 ) {
	// Cathode
	if ( k == 0 ) {
	  E_grid_r[j*NZ] = 0.;
	  E_grid_z[j*NZ] = 0.75*Phi[j*NZ] - Phi[j*NZ+1] + 0.25*Phi[j*NZ+2]; 
	}
	
	// Anode
	else if ( k == nz ) {
	  E_grid_r[j*NZ+nz] = 0.;
	  E_grid_z[j*NZ+nz] = -0.75*Phi[j*NZ+nz] + Phi[j*NZ+(nz-1)] - 0.25*Phi[j*NZ+(nz-2)]; 
	}
	
	// Symmetry axis
	else if ( (j == 0) && (0 < k) && (k < nz) ) {
	  E_grid_r[k] = 0.; // symmetry condition
	  E_grid_z[k] = 0.25* (Phi[k-1] - Phi[k+1]);
	}
	
	// Infinity
	else if ( (j == nr) && (0 < k) && (k < nz) ) {
	  E_grid_r[nr*NZ+k] = -0.75*Phi[nr*NZ+k] + Phi[(nr-1)*NZ+k] - 0.25*Phi[(nr-2)*NZ+k];
	  E_grid_z[nr*NZ+k] = 0.; 
	}
	
	// Inside
	else {
	  E_grid_r[j*NZ+k] = 0.25* (Phi[(j-1)*NZ+k] - Phi[(j+1)*NZ+k]);
	  E_grid_z[j*NZ+k] = 0.25* (Phi[j*NZ+(k-1)] - Phi[j*NZ+(k+1)]);
	}
      }
      
      //-----------------------------------------------------------------//
      else if ( BC == 1 || BC == 4 ) {
	// Cathode
	if ( k == 0 ) {
	  E_grid_r[j*NZ] = 0.;
	  E_grid_z[j*NZ] = 0.75*Phi[j*NZ] - Phi[j*NZ+1] + 0.25*Phi[j*NZ+2]; 
	}
	
	// Anode
	else if ( k == nz ) {
	  E_grid_r[j*NZ+nz] = 0.;
	  E_grid_z[j*NZ+nz] = -0.75*Phi[j*NZ+nz] + Phi[j*NZ+(nz-1)] - 0.25*Phi[j*NZ+(nz-2)]; 
	}
	
	// Symmetry axis
	else if ( (j == 0) && (0 < k) && (k < nz) ) {
	  E_grid_r[k] = 0.; // symmetry condition
	  E_grid_z[k] = 0.25* (Phi[k-1] - Phi[k+1]);
	}
	
	// Infinity
	else if ( (j == nr) && (0 < k) && (k < nz) ) {
	  E_grid_r[nr*NZ+k] = 0.;
	  E_grid_z[nr*NZ+k] = 0.25* (Phi[nr*NZ+(k-1)] - Phi[nr*NZ+(k+1)]); 
	}
	
	// Inside
	else {
	  E_grid_r[j*NZ+k] = 0.25* (Phi[(j-1)*NZ+k] - Phi[(j+1)*NZ+k]);
	  E_grid_z[j*NZ+k] = 0.25* (Phi[j*NZ+(k-1)] - Phi[j*NZ+(k+1)]);
	}
      }
      
      //-----------------------------------------------------------------//
      else if ( BC == 2 ) {
	// Cathode
	if ( (k == 0)  && (j > 0) && (j < nr) ) {
	  E_grid_r[j*NZ] = 0.25* (Phi[(j-1)*NZ] - Phi[(j+1)*NZ]);
	  E_grid_z[j*NZ] = 0.; 
	}
	
	// Anode
	else if ( (k == nz) && (j > 0) && (j < nr ) ) {
	  E_grid_r[j*NZ+nz] = 0.25* (Phi[(j-1)*NZ+nz] - Phi[(j+1)*NZ+nz]);
	  E_grid_z[j*NZ+nz] = 0.; 
	}
	
	// Rmin
	else if ( (j == 0) ) {
	  E_grid_r[k] = 0.75*Phi[k] - Phi[NZ+k] + 0.25*Phi[2*NZ+k]; 
	  E_grid_z[k] = 0.;
	}
	
	// Rmax
	else if ( (j == nr) ) {
	  E_grid_r[nr*NZ+k] = -0.75*Phi[nr*NZ+k] + Phi[(nr-1)*NZ+k] - 0.25*Phi[(nr-2)*NZ+k];
	  E_grid_z[nr*NZ+k] = 0.; 
	}
	
	// Inside
	else {
	  E_grid_r[j*NZ+k] = 0.25* (Phi[(j-1)*NZ+k] - Phi[(j+1)*NZ+k]);
	  E_grid_z[j*NZ+k] = 0.25* (Phi[j*NZ+(k-1)] - Phi[j*NZ+(k+1)]);
	}
      }
      
      //-----------------------------------------------------------------//
      else if ( BC == 3 ) {
	// Cathode
	if ( (k == 0)  && (j > 0) && (j < nr) ) {
	  //E_grid_r[j*NZ] = 0.25* (Phi[(j-1)*NZ] - Phi[(j+1)*NZ]);
	  E_grid_r[j*NZ] = 0.125* (Phi[(j-1)*NZ] + Phi[(j-1)*NZ+nz] 
				   - Phi[(j+1)*NZ] - Phi[(j+1)*NZ+nz]);
	  E_grid_z[j*NZ] = 0.25* (Phi[j*NZ+nz-1] - Phi[j*NZ+1]); 
	}
	
	// Anode
	else if ( (k == nz) && (j > 0) && (j < nr ) ) {
	  //E_grid_r[j*NZ+nz] = 0.25* (Phi[(j-1)*NZ+nz] - Phi[(j+1)*NZ+nz]);
	  E_grid_r[j*NZ+nz] = 0.125* (Phi[(j-1)*NZ] + Phi[(j-1)*NZ+nz] 
				      - Phi[(j+1)*NZ] - Phi[(j+1)*NZ+nz]);
	  E_grid_z[j*NZ+nz] = 0.25* (Phi[j*NZ+nz-1] - Phi[j*NZ+1]); 
	}
	
	// Symmetry axis
	else if ( (j == 0) ) {
	  E_grid_r[k] = 0.; 
	  
	  if ( k == 0 )
	    E_grid_z[k] = 0.25* (Phi[nz-1] - Phi[k+1]);
	  else if ( k == nz )
	    E_grid_z[k] = 0.25* (Phi[k-1] - Phi[1]);
	  else
	    E_grid_z[k] = 0.25* (Phi[k-1] - Phi[k+1]);
	}
	
	// Infinity
	else if ( (j == nr) ) {
	  E_grid_r[nr*NZ+k] = 0.;
	  if ( k == 0 )
	    E_grid_z[nr*NZ+k] = 0.25* (Phi[nr*NZ+(nz-1)] - Phi[nr*NZ+(k+1)]); 
	  else if ( k == nz )
	    E_grid_z[nr*NZ+k] = 0.25* (Phi[nr*NZ+(k-1)] - Phi[nr*NZ+1]); 
	  else
	    E_grid_z[nr*NZ+k] = 0.25* (Phi[nr*NZ+(k-1)] - Phi[nr*NZ+(k+1)]); 
	}
	
	// Inside
	else {
	  E_grid_r[j*NZ+k] = 0.25* (Phi[(j-1)*NZ+k] - Phi[(j+1)*NZ+k]);
	  E_grid_z[j*NZ+k] = 0.25* (Phi[j*NZ+(k-1)] - Phi[j*NZ+(k+1)]);
	}
      }

    }
  }
  
  // Integrate field for ions over dt_ion steps!
  for (int j=0; j<NR; j++) {
    for (int k=0; k<NZ; k++) {
      E_ion_z[j*NZ+k] += E_grid_z[j*NZ+k];
      E_ion_r[j*NZ+k] += E_grid_r[j*NZ+k];
    }
  }
}
