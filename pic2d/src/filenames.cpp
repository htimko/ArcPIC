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

  filenames.cpp:
  Creates filenames for data outputting 

***********************************************************************/

#include "pic.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define XTRN extern
#include "outp.h"
#undef XTRN



void  file_names_2D( const int nsteps ) {
  
  bool error = false;

  if ( snprintf(fn_e, LEN_FILENAME, "out/ne%08i.dat",   nsteps) >= LEN_FILENAME ) error=true ;
  if ( snprintf(fn_i, LEN_FILENAME, "out/Cup%08i.dat",  nsteps) >= LEN_FILENAME ) error=true ;
  if ( snprintf(fn_n, LEN_FILENAME, "out/Cu%08i.dat",   nsteps) >= LEN_FILENAME ) error=true ;
  
  if ( snprintf(fphi, LEN_FILENAME, "out/phi%08i.dat",  nsteps) >= LEN_FILENAME ) error=true ;
  
  if ( snprintf(fEz, LEN_FILENAME, "out/Ez%08i.dat",    nsteps) >= LEN_FILENAME ) error=true ;
  if ( snprintf(fEr, LEN_FILENAME, "out/Er%08i.dat",    nsteps) >= LEN_FILENAME ) error=true ;
  
  if ( snprintf(fv_ez, LEN_FILENAME, "out/uez%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;
  if ( snprintf(fv_er, LEN_FILENAME, "out/uer%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;
  if ( snprintf(fv_et, LEN_FILENAME, "out/uet%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;

  if ( snprintf(fT_ez, LEN_FILENAME, "out/Tez%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;
  if ( snprintf(fT_er, LEN_FILENAME, "out/Ter%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;
  if ( snprintf(fT_et, LEN_FILENAME, "out/Tet%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;

  if ( snprintf(fv_iz, LEN_FILENAME, "out/uiz%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;
  if ( snprintf(fv_ir, LEN_FILENAME, "out/uir%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;
  if ( snprintf(fv_it, LEN_FILENAME, "out/uit%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;
  
  if ( snprintf(fT_iz, LEN_FILENAME, "out/Tiz%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;
  if ( snprintf(fT_ir, LEN_FILENAME, "out/Tir%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;
  if ( snprintf(fT_it, LEN_FILENAME, "out/Tit%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;
  
  if ( snprintf(fr_e, LEN_FILENAME, "out/re%08i.dat",   nsteps) >= LEN_FILENAME ) error=true ;
  if ( snprintf(fr_i, LEN_FILENAME, "out/rCup%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;
  if ( snprintf(fr_n, LEN_FILENAME, "out/rCu%08i.dat",  nsteps) >= LEN_FILENAME ) error=true ;
  
  if ( snprintf(fvdf_ez, LEN_FILENAME, "out/vdfez%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;
  if ( snprintf(fvdf_er, LEN_FILENAME, "out/vdfer%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;
  if ( snprintf(fvdf_eabs, LEN_FILENAME, "out/vdfeabs%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;

  if ( snprintf(fvdf_iz, LEN_FILENAME, "out/vdfCupz%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;
  if ( snprintf(fvdf_ir, LEN_FILENAME, "out/vdfCupr%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;
  if ( snprintf(fvdf_iabs, LEN_FILENAME, "out/vdfCupabs%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;
  
  if ( snprintf(fvdf_nz, LEN_FILENAME, "out/vdfCuz%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;
  if ( snprintf(fvdf_nr, LEN_FILENAME, "out/vdfCur%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;
  if ( snprintf(fvdf_nabs, LEN_FILENAME, "out/vdfCuabs%08i.dat", nsteps) >= LEN_FILENAME ) error=true ;

  //Also have file from Arcbounds::writefile_removedParticles()

  if (error) {
    printf("\n\n");
    printf("Error detected in file_names_2D()\n");
    printf("Generated filenames:\n");
    printf("\t fn_e           : '%s'\n", fn_e);
    printf("\t fn_i           : '%s'\n", fn_i);
    printf("\t fn_n           : '%s'\n", fn_n);

    printf("\t fphi           : '%s'\n", fphi);

    printf("\t fv_ez          : '%s'\n", fv_ez);
    printf("\t fv_er          : '%s'\n", fv_er);
    printf("\t fv_et          : '%s'\n", fv_et);

    printf("\t fT_ez          : '%s'\n", fT_ez);
    printf("\t fT_er          : '%s'\n", fT_er);
    printf("\t fT_et          : '%s'\n", fT_et);

    printf("\t fv_iz          : '%s'\n", fv_iz);
    printf("\t fv_ir          : '%s'\n", fv_ir);
    printf("\t fv_it          : '%s'\n", fv_it);

    printf("\t fT_iz          : '%s'\n", fT_iz);
    printf("\t fT_ir          : '%s'\n", fT_ir);
    printf("\t fT_it          : '%s'\n", fT_it);

    printf("\t fr_e           : '%s'\n", fr_e);
    printf("\t fr_i           : '%s'\n", fr_i);
    printf("\t fr_n           : '%s'\n", fr_n);

    printf("\t fvdf_ez        : '%s'\n", fvdf_ez);
    printf("\t fvdf_er        : '%s'\n", fvdf_er);
    printf("\t fvdf_eabs      : '%s'\n", fvdf_eabs);

    printf("\t fvdf_iz        : '%s'\n", fvdf_iz);
    printf("\t fvdf_ir        : '%s'\n", fvdf_ir);
    printf("\t fvdf_iabs      : '%s'\n", fvdf_iabs);

    printf("\t fvdf_nz        : '%s'\n", fvdf_nz);
    printf("\t fvdf_nr        : '%s'\n", fvdf_nr);
    printf("\t fvdf_nabs      : '%s'\n", fvdf_nabs);
    
    exit(1);
  }

}
