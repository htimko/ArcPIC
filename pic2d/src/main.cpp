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
  
  main.cpp:
  Main program control and time stepping loop.
  Writing of LOCKFILE, out/timeIndex.dat, mainstats.dat

***********************************************************************/

#include  <stdio.h>
#include  <math.h>
#include  <sys/stat.h>
#include  <time.h>
#include  <string.h>
#include  <stdlib.h>

#include  <slu_ddefs.h>

#include <omp.h>

#define   XTRN
#include  "pic.h"
#include  "var.h"
#include  "dim.h"
#include  "arrays1.h"
#include  "outp.h"
#include  "mydef.h"

#include  "random.h"

#include  "init.h"
#include  "my_mem.h"
#include  "my_time.h"
#include  "phi.h"
#include  "efield.h"
#include  "push.h"
#include  "density.h"
#include  "e_ion.h"
#include  "moms.h"
#include  "aver.h"
#include  "filenames.h"
#include  "outputz.h"
#include  "engy.h"
#include  "colls.h"
#include  "order.h"
#include  "input.h"
#include  "vdf.h"
#include  "backup.h"
#include  "print_par.h"
#include  "checkbounds.h"
#include  "weightPotential.h"


void print_time( double t );

int main () {

  //Safety: avoid unintended double-starts
  struct stat st;
  if (not stat("LOCKFILE", &st) ) {
    printf("main(): LOCKFILE exists, refusing to start.\n");
    exit(1);
  }
  FILE* lockfile = fopen("LOCKFILE", "w");
  fprintf(lockfile, "This is a lockfile - ArcPic will refuse to start as long as it's present.\n");
  fprintf(lockfile, "Goal: Avoid unintensional restarts which corrupts older data.\n");
  fclose(lockfile);

  // Diagnostics
  int check_stab = 0; // PIC stability ok: 0=yes, 1=no
  int check_dens = 0; // Density well resolved: 0=yes, 1=no
  double f_dens;      // Factor for density check
  if ( snprintf( ferr, LEN_FILENAME, "%s", "error.log") >= LEN_FILENAME ) {
    printf("Error detected when making filename for ferr, generated filename = '%s'\n", ferr);
    exit(1);
  }
  
  // SuperLU parameters
  SuperMatrix L_slu, U_slu;
  int *perm_c_slu;
  int *perm_r_slu;
  double *rhs_slu;
  
  //Time index file
  // (written every output step, used to correlate step number with simulation time
  // for the sake of plotting etc)
  FILE* timeIndex;
  //Main loop stats file (written every time step,
  // contains vital overall statistics such as particle counts)
  FILE* mainStats;
  
  // Initialise variables 
  double dvt;
  int NG = NGR*NGZ;
  n_aver = 0, n_aver_ion = 0, n_aver_diagn = 0;
  
  //TODO: This should be backed up
  double induced_cathode_prev = 0.0;
  
  // Collision sanity checks
  Vec3d mcheck = {0.0,0.0,0.0};        
  double echeck = 0.0;
  
  //Read input.txt
  input();
  //Sanity check of dim.h
  if (NGR != nr+1 || NGZ != nz+1) {
    printf("ERROR!! (NGR,NGZ) = (%i,%i), while (nr,nz)=(%i,%i)! (should be one smaller)\n", NGR, NGZ,nr,nz);
    printf("Aborting!\n");
    exit(1);
  }
  //Sanity check of injection steps
  if (n2inj_step % dt_ion != 0 or n2inj_step == 0) {
    printf("ERROR!! n2inj_step=%i doesn't occur on dt_ion=%i\n", n2inj_step, dt_ion);
    printf("Aborting!\n");
    exit(1);
  }
  if (i2inj_step % dt_ion != 0 or i2inj_step == 0) {
    printf("ERROR!! i2inj_step=%i doesn't occur on dt_ion=%i\n", i2inj_step, dt_ion);
    printf("Aborting!\n");
    exit(1);
  }
  
  //Initialize multithreading
  if ( numParaThreads > omp_get_num_procs() ||  numParaThreads > omp_get_thread_limit() ) {
    printf("Error in main.cpp during initialization: Can't have numParaThreads=%i > num_procs=%i or thread_limit=%i\n", numParaThreads, omp_get_num_procs(), omp_get_thread_limit());
    exit(1);
  }
  else if (numParaThreads < 1) {
    printf("Error in main.cpp during initialization: Can't have numParaThreads=%i < 1\n", numParaThreads);
    exit(1);
  }
  omp_set_num_threads(numParaThreads);
  
  //File names for continuation
  if ( snprintf( fbckp, LEN_FILENAME, "out/bckp_%03i.dat", CONTINUATION ) >= LEN_FILENAME ) {
    printf("Error detected when making filename for fbckp, generated filename = '%s'\n", fbckp);
    exit(1);
  }
  if ( CONTINUATION > 0 ) {	 
    if ( snprintf( frestart, LEN_FILENAME, "out/bckp_%03i.dat", CONTINUATION-1 ) >= LEN_FILENAME ) {
      printf("Error detected when making filename for frestart, generated filename = '%s'\n", frestart);
      exit(1);
    }
  }
  
  dr = dz;
  
  double time_start = omp_get_wtime(); //Real walltime [s], relative to some (thread-dependent) starting point

  time_t now_is = my_time();
  printf("*** Starting up 2D Arc-PIC code, date:  *** %s\n", ctime(&now_is) );  
  
  //********INITIALISING...********//
  Names[0] =  "H+ ";
  Names[1] =  "Cu+";
  Names[2] =  "Cu ";
  
  // STARTING NEW RUN
  if ( CONTINUATION == 0 ) {
    nstepsmin = 0;
    printf("-- Starting new run. (%d steps) -- \n", nstepsmin);
    
    allocate_arrays( nr, nz, &perm_c_slu, &perm_r_slu, &rhs_slu );
    
    //Initialize scaling
    calc_parameters_2D();
    
    M_ions[0] =  mi_over_me;                 //  H+
    M_ions[1] =  63.546/1.0073*mi_over_me-1; //  Cu+
    M_ions[2] =  63.546/1.0073*mi_over_me;   //  Cu
    
    for (int sort=0; sort < NSpecies; sort++) {
      cs_ions[sort] =  sqrt(M_ions[0]/M_ions[sort])*vi_0;
      vt_ions[sort] =  sqrt(M_ions[0]/M_ions[sort])*v_ti;
    }
    
    rng_initialize(RNGbaseSeed, numParaThreads);
    printf("\n");
    rng_printStatus();
    printf("\n");

    // External circuit elements
    circuit->init();    
    //Particle boundary conditions
    pbounds->init(nr, Zmin, Zmax, Rmax);
    //Initial particle distribution
    if (iParts != NULL) {
      iParts->init();
    }

    // Add initial conditions, outputting
    if ( BC == 2 || BC == 3 ) potential_factorise_BC23( nr, nz, NR, NZ, dr, dz, &L_slu, &U_slu, &perm_c_slu, &perm_r_slu );
    else                        potential_factorise_2D( nr, nz, NR, NZ, dr, dz, &L_slu, &U_slu, &perm_c_slu, &perm_r_slu );
    
    // Initialize densities, field, and sputtering arrays
    for (int i=0; i<NSpecies*NGR*NGZ; i++) {
      n_i[i]=0.;
    } 
    for (int i=0; i<Lastion*NGR*NGZ; i++) {
      E_ion_r[i]=0.;
      E_ion_z[i]=0.;
    } 
    for (int i=0; i<NGR*NGZ; i++) {
      n_e[i]=0.;
      E_grid_r[i]=0.;
      E_grid_z[i]=0.;
      phi[i]=0.;
    }
    
    // Initialise particles
    nr_e = 0, nr_i[0] = 0, nr_i[1] = 0, nr_i[2] = 0;
    
    //Initialize time index file
    timeIndex = fopen("out/timeIndex.dat", "w");
    fprintf(timeIndex, "# StepNum SimTime[ns]\n");
    fflush(timeIndex);
    
    mainStats = fopen("mainStats.dat", "w");
    fprintf(mainStats, "# StepNum SimTime[ns] StepTime[s] n_electrons n_ions n_neutrals [superparticles]\n");
    fflush(mainStats);

    // INITIAL PARTICLE DISTRIBUTION
    if (iParts != NULL) {
      iParts->inject_e(elec, nr_e);
      iParts->inject_n(ions + Lastion*NPART, nr_i[Lastion]);
      for (int sort=0; sort<NSpecies; sort++) {
	if (q_ions[sort] != 0.) {
	  iParts->inject_i(ions+sort*NPART, nr_i[sort], sort);
	}
      }
      if ( OUT_COORD == 0 ) {
	// Output coordinates: position and velocity
	file_names_2D( 0 ); 
	out_coords_2D( elec, nr_e, 1, Omega_pe, dz, fr_e );
	out_coords_2D( ions + NPART, nr_i[1], dt_ion, Omega_pe, dz, fr_i );
	out_coords_2D( ions + Lastion*NPART, nr_i[Lastion], dt_ion, Omega_pe, dz, fr_n );
      }
      
    }

    // INITIAL ENERGY
    kin_pot_en( elec, nr_e, ions + NPART, nr_i[1], ions + Lastion*NPART, nr_i[Lastion], 
		&En_e, En_i+1, En_i+Lastion, &En_p, &En_tot, 
		1./M_ions[1], 1./M_ions[2], phi, NR, NZ, Omega_pe, dz );
    printf( "...... Initial energy .................... \n" );  
    printf( "e- ......: ne = %10zu En_e = %-9.5f \n",  nr_e, En_e );
    for (int sort = 0; sort < NSpecies; sort++) {
      printf( "%s......: np = %10zu En_p = %-9.5f \n", Names[sort], nr_i[sort], En_i[sort]);
    }
    printf( "Pot. en..:               En_pot = %-9.5f \n",  En_p );
    printf( "...... Total: En_tot = %-9.6f ........ \n", En_tot);
    printf( "\n");
    
    // TO TEST
    density_2D( n_e,      nr, nz, NR, NZ, Vcell, elec,         qe, nr_e,    0, 1 );
    density_2D( n_i + NG, nr, nz, NR, NZ, Vcell, ions + NPART, qi, nr_i[1], 1, 2 ); 
    file_names_2D( 0 ); 
    out_dens_2D( n_e,      1, -1., nr, nz, NZ, Omega_pe, dr, dz, fn_e ); // NEW 25.8.2010
    out_dens_2D( n_i + NG, 1,  1., nr, nz, NZ, Omega_pe, dr, dz, fn_i ); // NEW 25.8.2010
    for (int i=0; i<NSpecies*NGR*NGZ; i++) n_i[i]=0.;
    for (int i=0; i<NGR*NGZ; i++)          n_e[i]=0.;
    
    // I. GET POTENTIAL
    file_names_2D( 0 );
    out_phi_2D(phi, 1, nr, nz, NZ, Omega_pe, dr, dz, fphi);
    if ( BC == 2 || BC == 3 ) potential_backsolve_BC23( nr, nz, NR, NZ, dz, circuit->getU0(), circuit->getUNz(), 
							phi, L_slu, U_slu, perm_c_slu, perm_r_slu, n_e, n_i + NG, &rhs_slu );
    else                        potential_backsolve_2D( nr, nz, NR, NZ, dz, circuit->getU0(), circuit->getUNz(), 
							phi, L_slu, U_slu, perm_c_slu, perm_r_slu, n_e, n_i + NG, &rhs_slu );
    
    printf("\n");
    
    file_names_2D( 1 );
    out_phi_2D(phi, 1, nr, nz, NZ, Omega_pe, dr, dz, fphi);
    
    // II. CALCULATE FIELD
    electric_field_2D( phi, E_grid_r, E_grid_z, E_ion_r, E_ion_z, nr, nz, NR, NZ );
      
    printf("\n");
  }
  else {
    //Temporary
#warning backup/restart disabled!
    printf("Restart facility not fully implemented!\n");
    exit(1);

    // CONTINUING PREVIOUS RUN
    allocate_arrays( nr, nz, &perm_c_slu, &perm_r_slu, &rhs_slu );
    
    read_data( frestart );
    re_init();
    
    timeIndex = fopen("out/timeIndex.dat", "a");
    mainStats = fopen("mainStats.dat", "a");
    
    // Check data read in:
    checkbounds_2D( elec, nr_e, Rmin, Rmax, Zmin, Zmax );
    for (int sort = 0; sort < NSpecies; sort++) {
      checkbounds_2D( ions + sort*NPART, nr_i[sort], Rmin, Rmax, Zmin, Zmax ); 
    }
    
    // Compute L and U matrices
    if ( BC == 2 || BC == 3 ) potential_factorise_BC23( nr, nz, NR, NZ, dr, dz, &L_slu, &U_slu, &perm_c_slu, &perm_r_slu );
    else                        potential_factorise_2D( nr, nz, NR, NZ, dr, dz, &L_slu, &U_slu, &perm_c_slu, &perm_r_slu );
    
    printf("-- Continuing from old run. (%d steps) -- \n", nstepsmin);
    // Continue with loop - electron push
  }
  
  init_reactions();
  print_parameters_2D();
  
  
  printf( "*** Beginning main loop *** \n" );
  print_time( omp_get_wtime()-time_start);
  printf("\n");   

  double step_time = omp_get_wtime(); // time for each step, written to mainStats.dat

  
  //********BEGIN LOOP********//
  for (nsteps=nstepsmin; nsteps<=nstepsmax; nsteps++) {
    
    // III. MOVE PARTICLES
    if ( MAGNETIC == 0 )push_magnetic_2D( elec, E_grid_r, E_grid_z, Bz_ext, Bt_ext, 1., nr_e, NZ );
    else                         push_2D( elec, E_grid_r, E_grid_z, nr_e,                     NZ );
    
    pbounds->remove_e( elec, nr_e);
    
    if ( nsteps/e2inj_step*e2inj_step  == nsteps ) {
      pbounds-> inject_e(elec, nr_e, E_grid_z);
    }
    
    // IV. CALCULATE DENSITIES
    density_2D( n_e, nr, nz, NR, NZ, Vcell, elec, qe, nr_e, 0, 1 );
    
    if ( nsteps/dt_ion*dt_ion == nsteps ) { //BEGIN IF ION STEP

      scalEion_2D( E_ion_r, E_ion_z, NR, NZ, dt_ion, M_ions );  
      
      for (unsigned int sort=0; sort<NSpecies; sort++ ) {
	if (q_ions[sort] != 0.) {
	  if ( MAGNETIC == 0 ) push_magnetic_2D( ions + sort*NPART, E_ion_r + sort*NG, E_ion_z + sort*NG, Bz_ext, Bt_ext, -1.*dt_ion/M_ions[sort], nr_i[sort], NZ );
	  else                          push_2D( ions + sort*NPART, E_ion_r + sort*NG, E_ion_z + sort*NG,                                          nr_i[sort], NZ );
	  
	  pbounds->remove_i( ions + sort*NPART, nr_i[sort], sort );	  
	}
	else {
	  push_neutrals_2D( ions + sort*NPART, nr_i[sort] );
	  pbounds->remove_n(ions+sort*NPART, nr_i[sort]);
	}
      }
      
      if ( nsteps/i2inj_step*i2inj_step  == nsteps ) {
	for (int sort=0; sort<NSpecies; sort++) {
	  if (q_ions[sort] != 0.) {
	    pbounds->inject_i(ions+sort*NPART, nr_i[sort], E_grid_z, sort);
	  }
	}
      }
      
      for ( int sort=0; sort<NSpecies; sort++ ) {
	if (q_ions[sort] != 0.) {
	  density_2D( n_i + sort*NG, nr, nz, NR, NZ, Vcell, ions + sort*NPART, qi, nr_i[sort], sort, Lastion );
	}
      }
      
      if ( nsteps/n2inj_step*n2inj_step  == nsteps ) {
	pbounds->inject_n(ions + Lastion*NPART, nr_i[Lastion], E_grid_z);
      }	

      if( nsteps >= nav_start ) {
	density_2D( n_i + Lastion*NG, nr, nz, NR, NZ, Vcell, ions + Lastion*NPART, qi, nr_i[Lastion], Lastion, NSpecies ); // neutrals: only for outputting

	aver_moments_2D( mom_ion + NG, n_aver_ion, ions + NPART, nr_i[1], nr, nz, NZ );
	aver_moments_2D( mom_ion + Lastion*NG, n_aver_ion, ions + Lastion*NPART, nr_i[Lastion], nr, nz, NZ ); // should be _SN()
	average_2D( n_i + NG, n_i_av + NG, NR, NZ, n_aver_ion ); // NEW: 25.8.2010
	average_2D( n_i + Lastion*NG, n_i_av + Lastion*NG, NR, NZ, n_aver_ion ); // NEW: 25.8.2010
	
	n_aver_ion++;
	//printf("Averaging ion moments \n");
      }
      
      initEion_2D( E_ion_r, E_ion_z, NR, NZ );

    } //END IF ION STEP

    // V. MCC COLLISIONS
    
    // Coulomb collisions
    if (ncoll_el>0 && nsteps/ncoll_el*ncoll_el == nsteps) {
      // e-e Coulomb collisions
      order_2D(elec, nr_e, e_order, nr, nz, NZ);
      coll_el_knm_2D(elec, e_order, nr, nz, NZ, 1., 0, &mcheck, &echeck, ncoll_el);
      
      // i-i Coulomb collisions
      order_2D(ions + NPART, nr_i[1], i_order + NG, nr, nz, NZ);
      coll_el_knm_2D(ions + NPART, i_order + NG, nr, nz, NZ, M_ions[1], 1, &mcheck, &echeck, ncoll_el);
    }  
    
    // Other collisions		
    if(ncoll_ion > 0 && nsteps/ncoll_ion*ncoll_ion == nsteps ) {
      
      order_2D(elec, nr_e, e_order, nr, nz, NZ);                 //TODO: This is often not needed!
      order_2D(ions + NPART, nr_i[1], i_order + NG, nr, nz, NZ); //TODO: This is often not needed!
      order_2D(ions + Lastion*NPART, nr_i[Lastion], i_order + Lastion*NG, nr, nz, NZ);   
      
      //elastic Cu+ Cu collisions
      coll_ion_neutral_noSP_2D( ions + Lastion*NPART, i_order + Lastion*NG, M_ions[Lastion], 
				ions + NPART, i_order + NG, M_ions[1], 
				nr, nz, NZ, React_Cup_Cu_el, &mcheck, &echeck );  
      
      //elastic Cu Cu collisions
      coll_n_n_2D( ions + Lastion*NPART, i_order + Lastion*NG, nr, nz, NZ, React_Cu_Cu, &mcheck, &echeck ); 
      
      //elastic el Cu collisions 
      coll_el_all_fake_2D( ions + Lastion*NPART, i_order + Lastion*NG, M_ions[Lastion], 
			   elec, e_order, nr, nz, NZ, React_Cu_el ); 
      
      // e + Cu = Cu+ + 2e 
      coll_el_neutrals_2D( ions + Lastion*NPART, nr_i+Lastion, i_order + Lastion*NG, M_ions[Lastion], 
			   elec, &nr_e, e_order, ions + NPART, nr_i+1, 
			   nr, nz, NZ, React_Cu_ion, &mcheck, &echeck );  
    }
    
    
    // VI. UPDATE POTENTIAL
    //printf("induced cathode charge elec = %g, ions = %g, deltaQ = %g ", inducedCharge_cathode(elec, nr_e, e_order), inducedCharge_cathode(ions+NPART, nr_i[1], i_order+NG), pbounds->getDeltaQ());
    double induced_cathode = inducedCharge_cathode(ions+NPART, nr_i[1], i_order+NG) - inducedCharge_cathode(elec, nr_e, e_order);
    double induced_cathode_delta = induced_cathode - induced_cathode_prev; //positive for negative particles leaving surface and positive coming towards it
    //printf("induced current = %g wall current = %g total current = %g\n", induced_cathode_delta, pbounds->getDeltaQ(), pbounds->getDeltaQ() + induced_cathode_delta);
    induced_cathode_prev = induced_cathode;
    //circuit->timestep(pbounds->getDeltaQ(), nsteps);
    circuit->timestep(pbounds->getDeltaQ()+induced_cathode_delta, nsteps);
    

    //Notify the pbounds if this is a output timestep
    pbounds->timestep(nsteps, nsteps >= nav_start && nsteps == nav_start+nav_time);
    

    
    if ( BC == 2 || BC == 3 ) potential_backsolve_BC23( nr, nz, NR, NZ, dz, circuit->getU0(), circuit->getUNz(), 
							phi, L_slu, U_slu, perm_c_slu, perm_r_slu, n_e, n_i + NG, &rhs_slu );
    else                        potential_backsolve_2D( nr, nz, NR, NZ, dz, circuit->getU0(), circuit->getUNz(),
							phi, L_slu, U_slu, perm_c_slu, perm_r_slu, n_e, n_i + NG, &rhs_slu );
    
    // VII. CALCULATE FIELD
    electric_field_2D( phi, E_grid_r, E_grid_z, E_ion_r, E_ion_z, nr, nz, NR, NZ );
    
    // VIII. DIAGNOSTICS AND OUTPUT
    step_time = omp_get_wtime()-step_time;
    fprintf(mainStats, "%08i %010e %f %zu %zu %zu \n",
	    nsteps, nsteps*Omega_pe*1e9/(56414.6*sqrt(n_ref)), step_time,
	    nr_e, nr_i[1], nr_i[Lastion]);
    fflush(mainStats);
    step_time = omp_get_wtime();

    if (nsteps >= diagn_start) {
      aver_diagn_2D( n_e, diagn_dens, elec, diagn_Te, diagn_ne, nr_e, n_aver_diagn, nr, nz, NR, NZ );
      n_aver_diagn++;
      
      if (nsteps == diagn_start+nav_time) {
	diagn_av_stability( diagn_dens, diagn_Te, diagn_ne, n_aver_diagn, -1., v_te, nr, nz, 
			    NR, NZ, Omega_pe, dr, dz, nsteps, &check_stab, ferr );
	if (check_stab == 1) {
	  //printf("*** ERROR: stability violated. Stopping *** \n");
	  //goto The_END; //STOP!!
	  printf("*** ERROR: stability violated. *** \n");
	  check_stab = 0;
	  fflush(stdout);
	}
	
	n_aver_diagn = 0;
	diagn_start += dt_diagn;
      }
    }

    if (nsteps >= nav_start) {
      if ( OUT_VDF == 0 ) {
	dvt = 6.66666667/v_te;
	vel_dst_along_2D( vdf_ez, vdf_er, vdf_eabs, elec,                 nr_e,          nr, nz, dvt, n_aver );   
	
	dvt = 6.66666667/cs_ions[1];
	vel_dst_along_2D( vdf_iz, vdf_ir, vdf_iabs, ions + NPART,         nr_i[1],       nr, nz, dvt, n_aver );
	
	dvt = 6.66666667/cs_ions[Lastion];
	vel_dst_along_2D( vdf_nz, vdf_nr, vdf_nabs, ions + Lastion*NPART, nr_i[Lastion], nr, nz, dvt, n_aver );  
      }
      
      aver_moments_2D( mom_el, n_aver, elec, nr_e, nr, nz, NZ );
      average_2D( phi,      phi_av, NR, NZ, n_aver );
      average_2D( n_e,      n_e_av, NR, NZ, n_aver ); // NEW: 25.8.2010
      average_2D( E_grid_z, E_av_z, NR, NZ, n_aver ); // E-field output
      average_2D( E_grid_r, E_av_r, NR, NZ, n_aver ); // E-field output
      
      n_aver++;
      //printf("Averaging electron moments \n");
      
      if (nsteps == nav_start) {
	printf("Start averaging .......................... \n");  
	printf("Omega_pe*t= %-6.0f, t=%f ns (%d steps) \n",
	       nsteps*Omega_pe, nsteps*Omega_pe*1e9/(56414.6*sqrt(n_ref)), nsteps);
      }
      if( nsteps == nav_start+nav_time ) {
	printf("*** Outputting now. (%d steps) *** \n", nsteps);
	print_time( omp_get_wtime()-time_start);
	printf("No. electrons: %zu, ions: %zu, neutrals: %zu \n", nr_e,nr_i[1],nr_i[Lastion]);
	printf("\n");
	
	fprintf(timeIndex, "%08i %010e\n", nsteps, nsteps*Omega_pe*1e9/(56414.6*sqrt(n_ref)) );
	fflush(timeIndex);
	
	file_names_2D( nsteps ); 
	out_dens_2D( n_e_av,              n_aver,    -1., nr, nz, NZ, Omega_pe, dr, dz, fn_e );
	out_dens_2D( n_i_av + NG,         n_aver_ion, 1., nr, nz, NZ, Omega_pe, dr, dz, fn_i );
	out_dens_2D( n_i_av + Lastion*NG, n_aver_ion, 1., nr, nz, NZ, Omega_pe, dr, dz, fn_n );
	
	out_phi_2D( phi_av, n_aver, nr, nz, NZ, Omega_pe, dr, dz, fphi );
	
	out_vels_2D( mom_el, nr, nz, NZ, cs*sqrt(M_ions[0]/M_ions[1]), dr, dz, fv_ez, fv_er, fv_et ); 
	out_temps_2D( mom_el, cs*sqrt(M_ions[0]/M_ions[1]), me_over_mi*M_ions[0]/M_ions[1], nr, nz, NZ, dr, dz, fT_ez, fT_er, fT_et);   
	
	out_vels_2D( mom_ion + NG, nr, nz, NZ, cs*dt_ion*sqrt(M_ions[0]/M_ions[1]), dr, dz, fv_iz, fv_ir, fv_it); 
	out_temps_2D( mom_ion + NG, cs*dt_ion*sqrt(M_ions[0]/M_ions[1]), 1., nr, nz, NZ, dr, dz, fT_iz, fT_ir, fT_it); 
	
	if ( OUT_VDF == 0 ) {
	  out_fv_along_2D( vdf_ez, vdf_er, vdf_eabs, nr, nz, fvdf_ez, fvdf_er, fvdf_eabs );
	  out_fv_along_2D( vdf_iz, vdf_ir, vdf_iabs, nr, nz, fvdf_iz, fvdf_ir, fvdf_iabs );  
	  out_fv_along_2D( vdf_nz, vdf_nr, vdf_nabs, nr, nz, fvdf_nz, fvdf_nr, fvdf_nabs );
	}
	
	if ( OUT_EFIELD == 0 ) {
	  out_efield_2D( E_av_z, E_av_r, n_aver, nr, nz, NZ, Omega_pe, dr, dz, fEz, fEr );
	}
	
	if ( OUT_COORD == 0 ) {
	  // Output coordinates: position and velocity
	  out_coords_2D( elec, nr_e, 1, Omega_pe, dz, fr_e );
	  out_coords_2D( ions + NPART, nr_i[1], dt_ion, Omega_pe, dz, fr_i );
	  out_coords_2D( ions + Lastion*NPART, nr_i[Lastion], dt_ion, Omega_pe, dz, fr_n );
	}
	
	kin_pot_en( elec, nr_e, ions + NPART, nr_i[1], ions + Lastion*NPART, nr_i[Lastion], 
		    &En_e, En_i+1, En_i+Lastion, &En_p, &En_tot, 
		    1./M_ions[1], 1./M_ions[2], phi, NR, NZ, Omega_pe, dz );
	printf( "...... Energy balance .................... \n" );  
	printf( "e- ......: ne = %10zu En_e = %-9.5f \n",  nr_e, En_e );
	for (int sort = 0; sort < NSpecies; sort++) {
	  printf( "%s......: np = %10zu En_p = %-9.5f \n", Names[sort], nr_i[sort], En_i[sort]);           
	}
	printf( "Pot. en..:               En_pot = %-9.5f \n",  En_p );
	printf( "...... Total: En_tot = %-9.6f ........ \n", En_tot );
	printf( "\n" );
	
	fflush( stdout );

	// Check densities
	f_dens = SQU(1./Omega_pe)/n_aver;
	for (int k=0; k<=nz; k++ ) {
	  if ( (n_i_av[NG+k] > 5./f_dens) || (n_e_av[k] < -5./f_dens) ) {
	    printf("ERROR: UNDERRESOLVED n_i_av = %.4e > %.4e or n_e_av = %.4e < %.4e for k = %d \n",
		   n_i_av[NG+k],5./f_dens,n_e_av[k],-5./f_dens,k);
	    fflush(stdout);
	    check_dens = 1;
	    goto The_END;
	  }
	}
	
	n_aver = 0;
	n_aver_ion = 0;
	
	// Backup data (if we are not under-resolved)
#warning backup/restart disabled!
	//if ( check_dens == 0 ) save_data( fbckp );
	
	nav_start += nav_dt;
      }
    }
    
  }
  //********END LOOP********//
  
  
  
  // The_END: save_data( fbckp );
 The_END: delete_arrays( &perm_c_slu, &perm_r_slu, &rhs_slu );

  delete circuit;
  delete pbounds;

  fclose(timeIndex);

  fclose(mainStats);

  printf("Done!.......................... \n");
  
  printf( "*** Total runtime *** \n" );
  print_time( omp_get_wtime()-time_start);
  printf("\n");   
  fflush(stdout);
  
  return 0;
}

void  print_time( double t ) {
  // t = time in seconds						
  double  tt, hour, min, sec;
  
  t   = (floor)(t+0.5);
  tt  = modf(t/3600., &hour);
  tt  = modf(tt*60.,  &min);
  sec = t-hour*3600.-min*60.;
  
  printf("Computation time: %6.0f s (%2.2i:%2.2i:%2.2i) \n", 
	 t, (int)hour, (int)min, (int)sec );
}
