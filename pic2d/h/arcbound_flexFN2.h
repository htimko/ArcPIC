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

  arcbounds_flexFN2.h:
  Particle boundary conditions,
  Fowler-Nordheim with beta and alpha
  setable within "rings" of arbirary rmin,rmax

***********************************************************************/

#ifndef ARCBOUND_FLEXFN2_H
#define ARCBOUND_FLEXFN2_H

#include "arcbounds.h"
#include <vector>

//Baseclass for FlexFN2
class FlexFN2 : public ArcRemover {
 public:
  //Inject electron superparticles into a ring between r1 and r2 [dz]
  // with a setable alpha and beta, using the same injection algorithm as FlexFN.
  void injectFNring(Particle pa[], size_t &np,
		    double alpha, double beta, double const Ez[],
		    double r1, double r2);

};

class FlexFN2_ring : public FlexFN2 {
 public:
  FlexFN2_ring(std::vector<char*>& option);

  virtual void init   (unsigned int nr, double zmin, double zmax, double rmax);
  virtual void re_init(unsigned int nr, double zmin, double zmax, double rmax);

  virtual void print_par() const;

  virtual void inject_e(Particle pa[], size_t &np, double const Ez[]);
  virtual void inject_i(Particle pa[], size_t &np, double const Ez[],
			unsigned int sort) {}; //NOP
  virtual void inject_n(Particle pa[], size_t &np, double const Ez[]) {}; //NOP

  //Save and restore backup data
  virtual void backup(FILE* file) {};
  virtual void restoreBackup(FILE* file) {};  

  virtual const char* getName() const { return "FlexFN2_ring"; };
 private:
  double ring_r1, ring_r2;
  double ring_alpha, ring_beta;
};

#endif
