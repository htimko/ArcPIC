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

  e_ion.cpp:
  Re-scales the electric field to an "ion electric field" 

***********************************************************************/

#include <stdio.h>

#include  "dim.h"



void scalEion_2D( double ei_r[], double ei_z[], int NR, int NZ, int dti, double M_ions[] )
{
  int i, j, k;
  double *eid0_r, *eid0_z, *eid_r, *eid_z; // Field of species 0 and other species

  // To test
  //printf("In e_ion.cpp, before scaling Eg_z[1]= %.3e \n", ei_z[1]);
  eid0_r = ei_r;
  eid0_z = ei_z;

  // Only positive ions => *(-1)*dt_ion*me/mion

  // First (nr+1)*(nz+1) elements of the E_ion array
  for (j=0; j<NR; j++) 
    {
      for (k=0; k<NZ; k++) 
	{
	  eid0_r[j*NZ+k] *= -dti;
	  eid0_z[j*NZ+k] *= -dti;	  
	}
    }

  // Rest of the array
  for (i=Lastion-1; i>0; i--)
    {
      eid_r = ei_r + i*NGR*NGZ;
      eid_z = ei_z + i*NGR*NGZ;
      for (j=0; j<NR; j++) 
	{
	  for (k=0; k<NZ; k++) 
	    {
	      eid_r[j*NZ+k] = eid0_r[j*NZ+k]/M_ions[i];
	      eid_z[j*NZ+k] = eid0_z[j*NZ+k]/M_ions[i];
	    }
	}
    }

  for (j=0; j<NR; j++) 
    {
      for (k=0; k<NZ; k++) 
	{
	  eid0_r[j*NZ+k] /= M_ions[0];
	  eid0_z[j*NZ+k] /= M_ions[0];	  
	}
    }



}




void initEion_2D( double ei_r[], double ei_z[], int NR, int NZ )
{
  int  j,k;

  for (j=0; j<NR; j++) 
    {
      for (k=0; k<NZ; k++) 
	{
	  ei_r[j*NZ+k] = 0;
	  ei_z[j*NZ+k] = 0;
	}
    }
}
