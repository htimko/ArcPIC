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

  arrays1.h:
  Define external arrays

***********************************************************************/

#ifndef ARRAYS1_H
#define ARRAYS1_H

XTRN Particle *elec;
XTRN Particle *ions;
XTRN Particle *temp;

XTRN const char *Names[NSpecies];   // Name of ions
XTRN size_t nr_i[NSpecies];         // Nr of ions
XTRN double q_ions[NSpecies];       // Charge of ions
XTRN double M_ions[NSpecies];       // Mass of ions
XTRN double cs_ions[NSpecies];      // Sound velocity
XTRN double vt_ions[NSpecies];      // Thermal velocity 

XTRN Moments mom_el[NGR*NGZ],
  mom_ion[NSpecies*NGR*NGZ];

XTRN double Vcell[NGR];     

// Densities
XTRN double n_e[NGR*NGZ], n_i[NSpecies*NGR*NGZ],  
  n_e_av[NGR*NGZ], n_i_av[NSpecies*NGR*NGZ]; 

XTRN double phi[NGR*NGZ], phi_av[NGR*NGZ];
XTRN double E_grid_r[NGR*NGZ], E_grid_z[NGR*NGZ];
XTRN double E_av_r[NGR*NGZ], E_av_z[NGR*NGZ];
XTRN double E_ion_r[Lastion*NGR*NGZ], E_ion_z[Lastion*NGR*NGZ];

// Stability diagnosis
XTRN double diagn_Te[NGR*NGZ], diagn_ve[NGR*NGZ], 
  diagn_ne[NGR*NGZ], diagn_dens[NGR*NGZ];

// Ordering arrays for collisions
XTRN size_t e_order[NGR*NGZ], i_order[NSpecies*NGR*NGZ]; 

// Energy outputting
XTRN double En_i[NSpecies];

// VDF arrays
XTRN double* vdf_ez;
XTRN double* vdf_er;
XTRN double* vdf_eabs;
XTRN double* vdf_iz;
XTRN double* vdf_ir;
XTRN double* vdf_iabs;
XTRN double* vdf_nz;
XTRN double* vdf_nr;
XTRN double* vdf_nabs;

#endif
