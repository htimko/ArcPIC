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

  engy.h:
  Header file for engy.cpp

***********************************************************************/

void energy_electrons_2D( Particle el[], size_t ne, double *We );

void energy_all_2D( Particle el[], size_t ne, double *We, double *Wf, 
		    const double E_r[], const double E_z[], int NR, int NZ, 
		    double Ndb, double omega_pe, double dz);
 
void energy_total_2D( Particle el[], size_t ne, Particle ion[], size_t ni, Particle neu[], size_t nn,
		      double *We, double *Wi, double *Wn, double *Wf, double *Wt,
		      double me_over_MCup, double me_over_MCu,
		      const double E_r[], const double E_z[], int NR, int NZ, 
		      double Vcell[], double Ndb, double omega_pe, double dz );

void kin_pot_en( Particle el[], size_t ne, Particle ion[], size_t ni, Particle neu[], size_t nn,
		 double *We, double *Wi, double *Wn, double *Wp, double *Wt,
		 double me_over_MCup, double me_over_MCu,
		 const double Phi[], int NR, int NZ, double omega_pe, double dz );
