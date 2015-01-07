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

  weightPotential.h:  
  Calculate the charge induced on the electrodes
  using the Shockley-Ramo theorem.

  Returns the induced charge as a positive number.

***********************************************************************/

#ifndef WEIGHTPOTENTIAL_H
#define WEIGHTPOTENTIAL_H

#include "pic.h"
#include "dim.h"

#include <cstdlib>

double inducedCharge_cathode(Particle pa[], size_t np, size_t ordcount[]);
double inducedCharge_anode(Particle pa[], size_t np, size_t ordcount[]);

#endif
