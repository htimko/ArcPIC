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

  moms.h:
  Header file for moms.cpp

***********************************************************************/

void aver_moments_2D( Moments mom[], int n_av,
		      Particle pa[], size_t np,
		      int nr, int nz, int NZ );

void aver_diagn_2D( const double dens[], double dens_av[],
		    Particle pa[], double temp_av[], double np_av[],
		    size_t np, int n_av, int nr, int nz, int NR, int NZ );

void aver_moments_SN_2D( Moments mom[], int n_av,
			 Particle pa[], size_t np,
			 int nr, int nz, int NZ );


