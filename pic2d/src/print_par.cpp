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

  print_par.cpp:
  Prints read-in parameters to output file

***********************************************************************/

#include  "pic.h"
#include  "dim.h"

#include  <time.h>
#include  <stdio.h>
#include  <math.h>

#define   XTRN  extern 

#include  "var.h"
#include  "outp.h"
#include  "arrays1.h"



void print_parameters_2D( void )
{
        
  printf( "Input parameters initialised: \n" );
  printf( "\n" );
  printf( "_______________________________________________________________\n" );
  printf( " \n");
  printf( "Scaling and main parameters: \n" );
  printf( " - n_ref:               %g (in 1/cm^3)\n", n_ref );
  printf( " - T_ref:               %g (in eV)\n", T_ref );
  printf( " - Ndb:                 %g\n", Ndb );
  printf( " - nr:                  %d\n", nr );
  printf( " - nz:                  %d\n", nz );
  printf( " - dz:                  %g (in L_Db)\n", dz );
  printf( " - dt:                  %g (in O_pe^-1)\n", Omega_pe );
  printf( " - Ti_over_Te:          %g\n", Ti_over_Te );
  printf( " - mi_over_me:          %g\n", mi_over_me );
  //printf( " - seed:                %llu ( Note: may be changed by initrand() )\n", seed );
  printf( " - RNGbaseSeed          %lu\n", RNGbaseSeed );
  printf( " - numParaThreads       %i\n", numParaThreads );

  printf( " \n" );
  printf( "Timesteps: \n" );
  printf( " - dt_ion:              %d (in dt)\n", dt_ion );
  printf( " - ncoll_el:            %d (in dt)\n", ncoll_el );
  printf( " - ncoll_ion:           %d (in dt)\n", ncoll_ion );
  printf( " - dt_diagn:            %d (in dt)\n", dt_diagn );
  printf( " - nstepsmin:           %d\n", nstepsmin );
  printf( " - nstepsmax:           %d\n", nstepsmax );
  printf( " - nav_start:           %d (in dt)\n", nav_start );
  printf( " - nav_time:            %d (in dt)\n", nav_time );
  printf( " - nav_dt:              %d (in dt)\n", nav_dt );

  printf( " \n" );
  printf( "Fields, particles and boundary conditions: \n" );
  printf( " - Bz_ext:              %g\n", Bz_ext );
  printf( " - Bt_ext:              %g\n", Bt_ext );
  printf( " - e2inj_step:          %d (in dt)\n", e2inj_step );
  printf( " - n2inj_step:          %d (in dt)\n", n2inj_step );
  printf( " - i2inj_step:          %d (in dt)\n", i2inj_step );

  printf("\n");
  printf("Particle boundary parameters:\n");
  printf(" - ArcName: %s\n", pbounds->getName());
  pbounds->print_par();

  printf("\n");
  printf("Circuit parameters:\n");
  printf(" - CircuitName: %s\n", circuit->getName());
  circuit->print_par();

  printf("\n");
  printf("Initial particle distribution:\n");
  if (iParts == NULL) {
    printf(" -- No initial particle distribution initialized --\n");
  }
  else {
    printf(" - Initial particle distribution name: %s\n", iParts->getName());
    iParts->print_par();
  }

  printf( " \n" );
  printf( "Velocities: \n" );
  printf( " - v_te:                %g (in dt/dz)\n", v_te );
  printf( " - v_ti:                %g (in dt_ion/dz)\n", v_ti );
  printf( " - cs:                  %g (in dt/dz)\n", cs );
  printf( " - vi_0:                %g (in dt_ion/dz)\n", vi_0 );

  printf( " \n" );
  printf( "Options (0=yes 1=no): \n" );
  printf( " - Output coordinates:  %d \n", OUT_COORD );
  printf( " - Output vdf:          %d \n", OUT_VDF );
  printf( " - Magnetic push:       %d \n", MAGNETIC );
  printf( " - Continuing old run:  %d \n", CONTINUATION );

  printf( "\n" );
  printf( " - Field boundary condition BC = %i \n", BC );
  if (BC == 0)
    printf( "   (Phi=0 at r=nr aka infinity)\n");
  else if (BC == 1)
    printf( "   (DPhi/Dr=0 in r=nr aka infty)\n");
  else if (BC == 2)
    printf( "   (Phi const at r-boundaries, DPhi/Dz=0 at electrodes)\n");
  else if (BC == 3)
    printf( "   (Periodic B.C.)\n");
  else if (BC == 4)
    printf( "   (DPhi/Dr=0 in r=nr aka infty, alternative implementation)\n");

  printf( "\n" );
  printf( "\n" );

  printf( "Array allocation parameters: \n" );
  printf( " NPART = %i  NGZ = %i  NGR = %i  \n", NPART, NGZ, NGR );
  printf( "\n" );
  printf( "_______________________________________________________________\n" );
  printf( " \n");
  fflush(stdout);

}

      
