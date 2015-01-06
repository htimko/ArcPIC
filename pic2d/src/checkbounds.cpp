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

  checkbounds.cpp:
  checking read-in data for 2D boundary conditions

***********************************************************************/

#include  "pic.h"

#include  <stdio.h>



void checkbounds_2D( Particle pa[], size_t np, double rmin, double rmax, double zmin, double zmax )
{


  for (size_t n=0; n<np; n++)
    {

      if ( pa[n].p.r >= rmax )
	{ 
	  printf("Erroneous data: n=%zu  r=%g  rmax=%g\n", n, pa[n].p.r, rmax);
	  pa[n].p.r = rmax - 1.e-10;
	}

      else if ( pa[n].p.r < rmin )
	{ 						
	  printf("Erroneous data: n=%zu  r=%g  rmin=%g\n", n, pa[n].p.r, rmin);
	  pa[n].p.r = rmin + 1.e-10;
	}


      if ( pa[n].p.z >= zmax )
	{ 
	  printf("Erroneous data: n=%zu  z=%g  zmax=%g\n", n, pa[n].p.z, zmax);
	  pa[n].p.z = zmax - 1.e-10;
	}

      else if ( pa[n].p.z < zmin )
	{ 						
	  printf("Erroneous data: n=%zu  z=%g  zmin=%g\n", n, pa[n].p.z, zmin);
	  pa[n].p.z = zmin + 1.e-10;
	}
       
    }

  fflush(stdout);

}
         
