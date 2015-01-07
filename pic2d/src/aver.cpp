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
  
  aver.cpp
  Averages any function f=f[j*(nz+1)+k]

***********************************************************************/



void average_2D( const double gv[], double gv_av[], int NR, int NZ, int n_aver ) {
  if (n_aver == 0) {
    for (int j=0; j<NR; j++) {
      for (int k=0; k<NZ; k++) {
	gv_av[j*NZ+k] = gv[j*NZ+k];
      }
    }
  }
  else {
    for (int j=0; j<NR; j++) {
      for (int k=0; k<NZ; k++) {
	gv_av[j*NZ+k] += gv[j*NZ+k];
      }
    }
  }
}
