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

  arcbounds_sumTipFNField.h:
  Particle boundary conditions,
  Fowler-Nordheim getting the field at r=0
  by summing over all electrons directly,
  not by using the grid

***********************************************************************/

#ifndef ARCBOUND_SUMTIPFNFIELD
#define ARCBOUND_SUMTIPFNFIELD

#include "arcbounds.h"
#include <vector>

class SumTipFNField : public ArcRemover {
 public:

  SumTipFNField(std::vector<char*>& options);
  virtual ~SumTipFNField();
  virtual void print_par() const;
  
  //Inject particles
  virtual void inject_e(Particle pa[], size_t& np, double const Ez[]);
  virtual void inject_i(Particle pa[], size_t &np, double const Ez[], unsigned int sort) {}; // NOP
  virtual void inject_n(Particle pa[], size_t &np, double const Ez[]) {}; // NOP

  
  //Save and restore backup data
  virtual void backup(FILE* file) {}; //NOP
  virtual void restoreBackup(FILE* file) {}; //NOP

  virtual const char* getName() const { return "SumTipFNField"; };

  //Write to the and arcbounds_original.dat output files
  virtual void writeFile_extras(unsigned int nstep);
  FILE* ofile_tip;
  double field_PIC;
  double field_direct;
  double field_direct_particlesOnly;
  size_t tip_emitted;
  
 private:
  double beta_tip;
  double Remission, Remission_theor;
};

#endif
