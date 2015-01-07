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

  weightPotential.cpp:  
  Calculate the charge induced on the electrodes
  using the Shockley-Ramo theorem

***********************************************************************/

#include "weightPotential.h"

#define   XTRN extern
#include "var.h"
#undef    XTRN

double inducedCharge_cathode(Particle pa[], size_t np, size_t ordcount[]) {
  //TODO: It may be more accurate to use ordcount to group particles with similar weight,
  //starting with those with small weight.
  
  double induced = 0.0;
  for (size_t i = 0; i < np; i++) {
    induced += 1.0 - pa[i].p.z/nz;
  }
  
  return induced;
}

double inducedCharge_anode(Particle pa[], size_t np, size_t ordcount[]) {
  //NOT YET IMPLEMENTED!
  return 0.0;
}
