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

  backup.cpp:
  Backs up & reads back data for re-start cases

***********************************************************************/

#include  "dim.h"
#include  "pic.h"
#include  <string.h>
#include  <stdio.h>

#define   XTRN  extern
#include  "arrays1.h"
#include  "var.h"
#include  "outp.h"
#include  "mydef.h"
#include  "circuit.h"
#undef XTRN

void save_data( const char *fdata ) {
  
  FILE *file;
  file = fopen(fdata,"w");
  
  Particle *ion_sort;
  
  // INPUT AND CALC_PARAMETERS  
  //fprintf(file,"%19.11e %19.11e %19.11e %19.11e %19.11e %19.11e\n", n_ref, T_ref, Ndb, dz, dr, Omega_pe);
  fprintf(file,"%19.11e %19.11e %19.11e %19.11e %19.11e\n", n_ref, T_ref, dz, dr, Omega_pe); // possibility to re-start with different Ndb
  fprintf(file,"%19.11e %19.11e %19.11e %19.11e\n", Rmin, Rmax, Zmin, Zmax);
  fprintf(file,"%11d %11d %11d %11d %11zu\n", nr, NR, nz, NZ, nr_e);
  fprintf(file,"%11d %11d %11d %11d %11d %11d %19.11e\n",
	  nsteps, nstepsmax, dt_ion, dt_diagn, ncoll_el, ncoll_ion, Ampl);
  fprintf(file,"%11d %11d %11d\n",
	  e2inj_step, n2inj_step, i2inj_step);
  fprintf(file,"%11d %11d %11d %11d\n", nav_dt, nav_start, nav_time, diagn_start);
  fprintf(file,"%11d %11d %11d %11d\n", OUT_COORD, OUT_EFIELD, OUT_VDF, MAGNETIC);
  
  fprintf(file,"%19.11e %19.11e\n", Bz_ext, Bt_ext);
  fprintf(file,"%19.11e %19.11e %19.11e %19.11e\n", 
	  Ti_over_Te, mi_over_me, me_over_mi, lambda_De);
  fprintf(file,"%19.11e %19.11e %19.11e %19.11e\n", 
	  v_te, v_ti, cs, vi_0);
  
  circuit->backup(file);
  
  // ION ARRAYS	
  for(int sort=0; sort < NSpecies; sort++) {
    fprintf(file,"%11zu %23.15e %23.15e %23.15e %23.15e\n",
	    nr_i[sort], q_ions[sort], M_ions[sort], cs_ions[sort], vt_ions[sort]);
  }
  
  // VCELL
  for (int i=0; i<NR; i++) {
    fprintf(file,"%19.11e\n",Vcell[i]);
  }
  
  // PARTICLES
  for(size_t i=0; i < nr_e; i++) {
    fprintf(file,"%23.15e %23.15e %23.15e %23.15e %23.15e %5d %3d\n",
	    elec[i].p.z, elec[i].p.r, elec[i].p.vz,
	    elec[i].p.vr, elec[i].p.vt, elec[i].p.m, -1);
  }
  
  for(int sort=0; sort < NSpecies; sort++) {
    ion_sort = ions + sort*NPART;
    for(size_t i=0; i < nr_i[sort]; i++) {
      fprintf(file,"%23.15e %23.15e %23.15e %23.15e %23.15e %5d %3d\n",
	      ion_sort[i].p.z, ion_sort[i].p.r, ion_sort[i].p.vz,
	      ion_sort[i].p.vr, ion_sort[i].p.vt, ion_sort[i].p.m, sort);
    }
  }
  
  // E_GRID TO START FROM PUSHER
  for(int i=0; i < NR*NZ; i++) {
    fprintf(file,"%23lg %23lg", E_grid_z[i], E_grid_r[i]);
  }
  
  fprintf(file,"\n");
  
  fclose(file);
  
  printf("*** Saved to file '%s'. (%d steps) *** \n", fdata, nsteps);
  fflush(stdout);
  
}
      

void read_data( const char *fdata ) { 

  FILE *file;
  file = fopen(fdata,"r");

  int dummy;
  Particle *ion_sort;
  
  // INPUT AND CALC_PARAMETERS
  //fscanf(file,"%19lg %19lg %19lg %19lg %19lg %19lg", &n_ref, &T_ref, &Ndb, &dz, &dr, &Omega_pe);
  fscanf(file,"%19lg %19lg %19lg %19lg %19lg", &n_ref, &T_ref, &dz, &dr, &Omega_pe); // possibility to re-start with different Ndb
  fscanf(file,"%19lg %19lg %19lg %19lg", &Rmin, &Rmax, &Zmin, &Zmax);
  fscanf(file,"%11d %11d %11d %11d %zu", &nr, &NR, &nz, &NZ, &nr_e);
  fscanf(file,"%11d %11d %11d %11d %11d %11d %19lg",
	 &nsteps, &nstepsmax, &dt_ion, &dt_diagn, &ncoll_el, &ncoll_ion, &Ampl);
  fscanf(file,"%11d %11d %11d", &e2inj_step, &n2inj_step, &i2inj_step);
  fscanf(file,"%11d %11d %11d %11d", &nav_dt, &nav_start, &nav_time, &diagn_start);
  fscanf(file,"%11d %11d %11d %11d", &OUT_COORD, &OUT_EFIELD, &OUT_VDF, &MAGNETIC);
  
  fscanf(file,"%19lg %19lg", &Bz_ext, &Bt_ext);
  fscanf(file,"%19lg, %19lg %19lg %19lg", 
	 &Ti_over_Te, &mi_over_me, &me_over_mi, &lambda_De);
  fscanf(file,"%19lg %19lg %19lg %19lg", 
	 &v_te, &v_ti, &cs, &vi_0);
  
  circuit->restoreBackup(file);
  
  // ION ARRAYS	
  for (int sort=0; sort < NSpecies; sort++) {
    fscanf(file,"%zu %23lg %23lg %23lg %23lg",
	   nr_i+sort, q_ions+sort, M_ions+sort, cs_ions+sort, vt_ions+sort);
  }
  
  // VCELL
  for (int i=0; i<NR; i++){
    fscanf(file,"%19lg",&(Vcell[i]));
  }
  
  // PARTICLES
  for(size_t i=0; i < nr_e; i++) {
    fscanf(file,"%23lg %23lg %23lg %23lg %23lg %5d %3d",
	   &(elec[i].p.z), &(elec[i].p.r), &(elec[i].p.vz),
	   &(elec[i].p.vr), &(elec[i].p.vt), &(elec[i].p.m), &dummy );
  }
  for(int sort=0; sort < NSpecies; sort++) {
    ion_sort = ions + sort*NPART;
    for(size_t i=0; i < nr_i[sort]; i++) {
      fscanf(file,"%23lg %23lg %23lg %23lg %23lg %5d %3d",
	     &(ion_sort[i].p.z), &(ion_sort[i].p.r), &(ion_sort[i].p.vz), 
	     &(ion_sort[i].p.vr), &(ion_sort[i].p.vt), &(ion_sort[i].p.m), &dummy );
    }
  }
  
  
  // E_GRID TO START FROM PUSHER
  for(int i=0; i < NR*NZ; i++) {
    fscanf(file,"%23lg %23lg", E_grid_z+i, E_grid_r+i);
  }
  
  fclose(file);
  
  printf("*** Re-starting from file %s. (%d steps) *** \n", fdata, nsteps);
  
}
