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

  dim.h:
  Defines memory allocation parameters

***********************************************************************/

#define  NGR      481       // nr+1
#define  NGZ      121       // nz+1

#define  NPART    200000000
#define  Lastion  2         // last ion with q != 0
#define  NSpecies 3         // non-electron species  

#define  Nvdst    401       // vdf resolution, choose odd nr
