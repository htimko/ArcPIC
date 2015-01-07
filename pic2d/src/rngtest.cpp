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

  rngtest.cpp:
  Test routine for checking that the RNG works properly

***********************************************************************/

#include <iostream>
#include <fstream>
using namespace std;

#include "random.h"

int main(int argc, char* argv[]) {
  unsigned long int seed = 1234;
  rng_initialize(seed, 1);
  
  int numCalls = 50000;
  
  //Test GaussRandom: mean = 0.5, sigma=1.5
  ofstream gaussFile;
  gaussFile.open("gausRandom.out");
  gaussFile.precision(10);
  for (int i = 0; i < numCalls; i++) {
    gaussFile << GausRandom(0.5,1.5) << " ";
  }
  gaussFile.close();
  
  //Wanted: Test the autocorrelation for RAND
}
