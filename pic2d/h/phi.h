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

  phi.h:
  Header file for phi.cpp

***********************************************************************/

void potential_factorise_2D( int nr, int nz, int NR, int NZ, double dr, double dz,
			     SuperMatrix* L, SuperMatrix* U, int** perm_c, int** perm_r );

void potential_factorise_BC23( int nr, int nz, int NR, int NZ, double dr, double dz,
			      SuperMatrix* L, SuperMatrix* U, int** perm_c, int** perm_r );


void potential_backsolve_2D( int nr, int nz, int NR, int NZ, double dz,
			     double const phi0, double const phiNz, double Phi[],
			     SuperMatrix L, SuperMatrix U, int* perm_c, int* perm_r,
			     double n_e[], double n_i[], double** rhs );

void potential_backsolve_BC23( int nr, int nz, int NR, int NZ, double dz,
			     double const phi0, double const phiNz, double Phi[],
			     SuperMatrix L, SuperMatrix U, int* perm_c, int* perm_r,
			      double n_e[], double n_i[], double** rhs );
