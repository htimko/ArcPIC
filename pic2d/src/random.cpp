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

  random.cpp:
  Interface to GSL random number generators

***********************************************************************/

#include "random.h"

#include <iostream>
#include <cstdlib>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

unsigned int numStreams_init = 0; //How many streams are initialized?
gsl_rng** random_streams = NULL;
unsigned long int* random_seeds = NULL; //The seeds used to initialize the RNG streams
unsigned long int random_baseSeed = 0;

void rng_initialize(unsigned long int baseSeed, unsigned int numStreams) {
  //Initialize a seed-generator
  if (baseSeed > 2147483648) {
    cout << "Expects: baseSeed < 2^31" << endl;
    exit(1);
  }
  gsl_rng* seedStream = gsl_rng_alloc( gsl_rng_rand );
  gsl_rng_set( seedStream, baseSeed );
  random_baseSeed = baseSeed;

  //Initialize the streams
  if (numStreams == 0) {
    cout << "Needs to initialize at least 1 RNG stream." << endl;
    exit(1);
  }
  numStreams_init = numStreams;
  random_streams = new gsl_rng*[numStreams];
  random_seeds = new unsigned long int[numStreams];
  for (unsigned int i = 0; i < numStreams; i++) {
    random_streams[i] = gsl_rng_alloc( RANDOM_RNG );
    random_seeds[i] = gsl_rng_get( seedStream );
    gsl_rng_set( random_streams[i], random_seeds[i] );
  }
  
  gsl_rng_free(seedStream); seedStream = NULL;
}
void rng_free() {
  if (numStreams_init == 0) {
    cout << "Can't call rng_free() when there are none to free" << endl;
    exit(1);
  }
  for (unsigned int i = 0; i < numStreams_init; i++) {
    gsl_rng_free( random_streams[i] );
    random_streams[i] = NULL;
  }
  delete [] random_streams; random_streams = NULL;
  delete[] random_seeds; random_seeds = NULL;
  numStreams_init = 0;
  random_baseSeed = 0;
}
unsigned int rng_getNumStreams() {
  return numStreams_init;
}
void rng_printStatus(bool hexDump) {
  cout << "RNG status:" << endl;
  cout << "Number of streams: " << numStreams_init << endl;
  cout << "Base seed: " << random_baseSeed << endl;
  cout << "Generator type: " << RANDOM_RNG->name << endl;
  cout << "Per-stream statuses:" << endl;
  for (unsigned int i = 0; i < numStreams_init; i++) {
    cout << "Stream #" << i << endl;
    cout << "Seed: " << random_seeds[i] << endl << flush;
    if (hexDump) {
      cout << "State: ";
      gsl_rng_print_state( random_streams[i] );
      cout << endl;
    }
  }
}

double Random()                      { return Random(0); };
double Random(unsigned int streamID) { return gsl_rng_uniform(random_streams[streamID]); };

double GausRandom(double mean, double sigma) {
  return GausRandom(mean, sigma, 0);
};
double GausRandom(double mean, double sigma, unsigned int streamID) {
  return gsl_ran_gaussian( random_streams[streamID], sigma ) + mean;
};
