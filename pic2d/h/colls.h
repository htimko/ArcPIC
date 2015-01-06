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

  colls.h:
  Header file for collisions.cpp

***********************************************************************/

#ifndef COLLS_H
#define COLLS_H

void coll_el_knm_2D( Particle  pa[], size_t ordcount[],
		     int nr,int nz,int NZ, 
		     double Mpa_over_me, int kind,
		     Vec3d *momcheck, double *engcheck, int ncoll);

void coll_ion_neutral_noSP_2D( Particle  neutrals[], size_t ordcount_ntrls[], double M_n,
			       Particle  ions[], size_t ordcount_ion[], double M_i,
			       int nr, int nz, int NZ, Reaction React,
			       Vec3d *momcheck, double *engcheck   );

#define COLL_N_N_2D_OMP_MINPARTICLES 1000 //Lower bound for where paralellization kicks in
void coll_n_n_2D( Particle neutrals[], size_t ordcount_ntrls[],                              
		  int nr, int nz, int NZ, Reaction React,
		  Vec3d *momcheck, double *engcheck );

void coll_el_all_fake_2D( Particle molecules[], size_t ordcount_m[], double M_m,   // molecules
			  Particle electrons[], size_t ordcount_el[],              // electrons
			  int nr, int nz, int NZ, Reaction React);

void coll_el_neutrals_2D( Particle  neutrals[], size_t *nn, size_t ordcount_ntrls[], double M_n,
			  Particle  electrons[], size_t *ne, size_t ordcount_el[],
			  Particle  ions[], size_t *ni, int nr, int nz, int NZ, Reaction React,
			  Vec3d *momcheck, double *engcheck );

#endif
