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

  var.h:
  Global variables

***********************************************************************/

#ifndef VAR_H
#define VAR_H

#include "circuit.h"
#include "arcbounds.h"
#include "initialParticles.h"

XTRN double Ampl,
  av_start, av_time,
  Bt_ext, Bz_ext,
  cs,
  dr, dz,
  dt_out,
  lambda_De,
  me_over_mi,
  mi_over_me,
  Ndb,
  n_ref, T_ref,       /*  Ref. dens. and temp. for rescaling  */
  Omega_pe,
  qe, qi,
  Ti_over_Te,
  v_te, v_ti,
  vi_0,
  Zmin, Zmax, Rmin, Rmax;

XTRN Circuit* circuit;
XTRN ArcBounds* pbounds;
XTRN InitialParticles* iParts;

XTRN int diagn_start,
  dt_diagn, 
  dt_ion,
  n_aver,
  n_aver_diagn,
  n_aver_ion,
  nav_start,
  nav_time,
  nav_dt,
  ncoll_el, ncoll_ion,
  nr, nz,
  NR, NZ,
  nsteps, nstepsmax, nstepsmin;

//Array indices/used-lengths
XTRN size_t nr_e;

XTRN int CONTINUATION,
  MAGNETIC,
  OUT_COORD,
  OUT_EFIELD,
  OUT_VDF,
  BC;

XTRN int e2inj_step, i2inj_step, n2inj_step;

XTRN double En_e, En_f, En_p, En_tot; //En_i, En_n, 

XTRN Reaction  React_Cu_el, React_Cu_ion, React_Cup_Cu_el, React_Cu_Cu;
XTRN double Ldb, N_sp;

//RNG and paralellization setup
XTRN unsigned long int RNGbaseSeed;
XTRN int numParaThreads;

#endif
