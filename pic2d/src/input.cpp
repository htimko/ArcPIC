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

  input.cpp:
  Reads in user-given parameters from input.txt

***********************************************************************/

#include  "pic.h"

#define XTRN extern
#include  "var.h"
#include  "circuit.h"
#include  "arcbounds.h"
#undef XTRN


#include "filenames.h"
#include "input.h"

#include  <stdio.h>

#include <cstdlib>
#include <sys/stat.h>
#include <iostream>
//#include <string.h>
#include <cstdlib>
#include <cstring>

using namespace std;

void input( void ) {

  FILE *in_file;

  //Check that file exists, else we get a segfault
  struct stat st;
  if ( stat("input.txt", &st) ) {
    printf("input(): Could not find file \"input.txt\". Aborting!\n");
    exit(1);
  }
  //Also check that output folder exist, create if neccessary
  if ( stat("out", &st) ) {
    printf("input(): Folder \"out\" did not exist. Creating it now...\n");
    if ( mkdir("out", S_IRWXU) ) {
      printf("input(): Could not create folder \"out\". Abort!\n");
      exit(1);
    }
  }
  
  in_file = fopen("input.txt","r");

  //SCALING PARAMETERS
  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%lg", &n_ref);
  
  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%lg", &T_ref);
  
  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%lg", &Ndb);
  
  
  fscanf(in_file,"%*[^:]%*[:]");
  fscanf(in_file,"%d %d", &nr, &nz);
  
  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%lg", &dz);
  
  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%lg", &Omega_pe);
  
  
  // TIMESTEPS
  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%d", &dt_ion);
  
  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%d", &ncoll_el);
  
  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%d", &ncoll_ion);
  
  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%d", &dt_diagn);
  
  
  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%d", &nstepsmax);
  
  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%lg", &dt_out);
  
  fscanf(in_file,"%*[^:]%*[:]");
  fscanf(in_file,"%lg", &av_start );
  
  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%lg", &av_time );
  
  
  // FIELDS, PARTILCES AND BOUNDARY CONDITIONS

  fscanf(in_file,"%*[^:]%*[:]");
  fscanf(in_file,"%lg", &Bz_ext);
  
  fscanf(in_file,"%*[^:]%*[:]");
  fscanf(in_file,"%lg", &Bt_ext);
  
  //Injection timesteps
  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%d", &e2inj_step);

  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%d", &n2inj_step);
  
  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%d", &i2inj_step);
  
  //New-style particle boundary conditions
  pbounds = ArcBounds::LoadArcBounds(in_file);

  // Circuit parameters
  circuit = Circuit::LoadCircuit(in_file); //Global var in var.h

  //Initial particles
  iParts = InitialParticles::LoadInitialParticles(in_file);

  // OPTIONS
  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%d", &OUT_COORD);
  
  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%d", &OUT_EFIELD);
  
  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%d", &OUT_VDF);
  
  
  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%d", &MAGNETIC);
  
  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%d", &CONTINUATION);
  
  
  // MISCELLANEOUS
  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%lg", &Ti_over_Te);
  
  fscanf(in_file,"%*[^:]%*[:]"); 
  fscanf(in_file,"%lg", &mi_over_me);
  
  fscanf(in_file,"%*[^:]%*[:]");
  fscanf(in_file,"%lu", &RNGbaseSeed);
  
  fscanf(in_file,"%*[^:]%*[:]");
  fscanf(in_file,"%i", &BC);
  
  fscanf(in_file,"%*[^:]%*[:]");
  fscanf(in_file,"%i", &numParaThreads);
  
  fclose (in_file);
  
}


char* readInputSection (FILE* in_file, vector<char*>& options_ret, bool acceptNone) {
  char wantName [NAME_MAXLEN];
  char* retName = NULL;
  vector<char*> options;
  
  fscanf(in_file, "%*[^:]%*[:]");
  fscanf(in_file, "%sNAME_MAXLEN", wantName);
  //Move the file pointer past the '\n'
  // foo should be much longer that NAME_MAXLEN (foo also used below)
  char foo[NAME_MAXLEN];
  memset(foo, '\0', NAME_MAXLEN);
  fgets(foo,NAME_MAXLEN,in_file);
  if (foo[NAME_MAXLEN-1] != '\0') {
    cout << "Error in readInputSection(): Ran out of space in foo (pos1)." 
	 << " Check whitespace in input.txt!" << endl;
    exit(1);
  }

  cout << "wantName='" << wantName << "'\n";

  //Loop over config sections
  const unsigned int safetyMax = 100;
  for (unsigned int safety = 0; safety < safetyMax; safety++) {
    //Loop sanity check
    if ( safety >= safetyMax ) {
      cout << "Error in readInputSection: Infinite loop in parser." << endl;
      exit(1);
    }
    
    //Find the right section
    char thisName [NAME_MAXLEN];
    memset(foo, '\0', NAME_MAXLEN);
    fgets(foo, NAME_MAXLEN,in_file);
    if (foo[NAME_MAXLEN-1] != '\0') {
      cout << "Error in readInputSection(): Ran out of space in foo (pos2)."
	   << " Check input.txt!" << endl;
      exit(1);
    }
    sscanf(foo, "*** %sNAME_MAXLEN", thisName);

    cout << "thisName='" << thisName << "'\n";    

    if ( ! strncmp(thisName, "END", NAME_MAXLEN) ) break;
    
    //Loop over options for thisName
    while(true) {
      //Get the next line
      char* line = new char[LINE_MAXLEN];
      memset(line, '\0', LINE_MAXLEN);
      fgets(line, LINE_MAXLEN, in_file);
      if ( line[LINE_MAXLEN-1] != '\0') {
	cout << "Error in readInputSection(): Ran out of space in line. "
	     << " Check input.txt!" << endl;
	exit(1);
      }
      //printf("Got line '%s'\n", line); //DEBUG

      //Check for "///"
      if ( ! strncmp(line, "///", 3) ) break;

      //Collect the line
      options.push_back(line);
    }
    
    //Got the right section!
    if ( ! strncmp( thisName, wantName, NAME_MAXLEN ) ) {
      cout << "In readInputSection(), got name = '" << wantName << "'" << endl;
      cout << "Got options:" << endl;
      for (size_t i = 0; i < options.size(); i++) {
	cout << "\t[" << i << "]:" << options[i]; //endl not needed, already \n at end of lines
      }
      
      if (retName != NULL) {
	cout << "Error in readInputSection: Found the wanted name twice in input.txt!" << endl;
	exit(1);
      }
      
      //Don't return immediatly, need to parse the rest of the sections
      retName = new char[NAME_MAXLEN];
      strncpy(retName, wantName, NAME_MAXLEN);
      for (size_t i = 0; i < options.size(); i++) {
	char* line_ret = new char[LINE_MAXLEN];
	strncpy(line_ret, options[i], LINE_MAXLEN);
	options_ret.push_back(line_ret);
      }
    }
    
    //Wrong section, clear the options array before next iteration
    for (size_t i = 0; i < options.size(); i++) {
      delete[] options[i];
    }
    options.clear();
  }
  
  if (retName != NULL) {
    //Implicitly also return options_ret
    return retName;
  }
  else if ( ! strncmp( wantName, "None", NAME_MAXLEN ) and acceptNone ) {
    retName = new char[NAME_MAXLEN];
    strncpy(retName, wantName, NAME_MAXLEN);
    return retName;
  }
  //Implicit else
  cout << "Error in readInputSection(): Couldn't find section corresponding to "
       << " '" << wantName << "'" << endl;
  exit(1);
  return NULL; //avoid compiler warning
}
