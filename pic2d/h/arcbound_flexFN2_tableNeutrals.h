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

  arcbounds_flexFN2_tableNeutrals.h:
  Particle boundary conditions,
  Neutral sputtering from tabulated distrubutions,
  and FlexFN2 electron injection

***********************************************************************/

#ifndef ARCBOUND_FLEXFN2_TABLENEUTRALS_H
#define ARCBOUND_FLEXFN2_TABLENEUTRALS_H

#include "arcbound_flexFN2.h"
#include <vector>

class FlexFN2_tableNeutrals : public FlexFN2 {
 public:
  FlexFN2_tableNeutrals(std::vector<char*>& option);

  virtual void init   (unsigned int nr, double zmin, double zmax, double rmax);
  virtual void re_init(unsigned int nr, double zmin, double zmax, double rmax);

  virtual void print_par() const;

  virtual void inject_e(Particle pa[], size_t &np, double const Ez[]);
  virtual void inject_i(Particle pa[], size_t &np, double const Ez[],
			unsigned int sort) {}; //NOP
  virtual void inject_n(Particle pa[], size_t &np, double const Ez[]);

  //Save and restore backup data
  virtual void backup(FILE* file) {};
  virtual void restoreBackup(FILE* file) {};  

  virtual const char* getName() const { return "FlexFN2_tableNeutrals"; };
 private:
  double emitter_radius;
  double emitter_alpha;
  double emitter_beta;

  double flat_alpha;
  double flat_beta;
};

#endif
