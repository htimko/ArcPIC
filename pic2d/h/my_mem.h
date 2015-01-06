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

  my_mem.h:
  Header file for my_mem.cpp

***********************************************************************/

void allocate_arrays( int nr, int nz, int** perm_c, int** perm_r, double** rhs );

void delete_arrays( int** perm_c, int** perm_r, double** rhs );
 
