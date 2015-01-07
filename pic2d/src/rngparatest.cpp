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

  rngparatest.cpp:
  Test RNG parallelization

***********************************************************************/


#include <iostream>
#include <cstdlib>

#include <omp.h>

#include "random.h"

using namespace std;

int main(int argc, char* argv[]) {
  if (argc != 3) {
    cout << "Usage: ./rngparatest #threads #iterations/thread" << endl;
    exit(1);
  }
  int num_threads = atoi(argv[1]);
  int num_iters = atoi(argv[2]);
  
  cout << "User requested:" << endl;
  cout << "num_threads = " << num_threads << endl;
  cout << "num_iters = " << num_iters << endl;
  cout << endl;

  if (num_threads > omp_get_thread_limit()) {
    num_threads = omp_get_thread_limit();
    cout << "Limited by thread limit." << endl;
  }
  omp_set_num_threads(num_threads);

  cout << "Number of threads :" << num_threads << endl;
  cout << "Number of CPUs    :" << omp_get_num_procs() << endl;
  cout << "Thread limit      :" << omp_get_thread_limit() << endl;

  //rng_printStatus();
  rng_initialize(1234, num_threads);
  rng_printStatus(true);
  cout << endl;
  
  cout << "Calculating " << num_iters << " iterations..."<< endl;

#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    for (int i = 0; i < num_iters; i++) {
      cout << tid << " " << i << " " << Random(tid) << endl << flush;
      //Random(tid);
    }
  }
  cout << " done." << endl << endl;

  rng_printStatus(true);

}
