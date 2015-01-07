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

  order.h:
  Header file for order.cpp

***********************************************************************/

//More memory-efficient ordering routine.
// Input:
// - pa: Particles that should be ordered.
//   May be provided in any order
// - np: Number of particles contained in pa (indices 0 -- np-1 valid)
// - nr,nz: Number of cells in r,z direction
// - NZ: Number of nodes in z direction = (nz+1)
// Output:
// - ordcount: Starting index of 
void order_2D(Particle pa[], size_t np, size_t ordcount[],
	      int nr, int nz, int NZ);
