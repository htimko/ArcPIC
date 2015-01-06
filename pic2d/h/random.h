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

  random.h:
  Interface to GSL random number generators,
  allowing for paralell streams

***********************************************************************/

#ifndef RANDOM_H
#define RANDOM_H

//Which of GSL's RNGs to use?
#define RANDOM_RNG gsl_rng_ranlxs2
//#define RANDOM_RNG gsl_rng_taus2 //Faster but not as well proven

#define RAND Random()

void rng_initialize(unsigned long int baseSeed, unsigned int numStreams);
void rng_free();
unsigned int rng_getNumStreams();
void rng_printStatus(bool hexDump=false);

//Get a uniform random number
double Random();
double Random(unsigned int streamID);

//Other distributions
double GausRandom(double mean, double sigma);
double GausRandom(double mean, double sigma, unsigned int streamID);

#endif
