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

  arcbounds_flexFN.h:
  Particle boundary conditions,
  Fowler-Nordheim with beta and alpha functions of r

***********************************************************************/

#ifndef ARCBOUND_FLEXFN_H
#define ARCBOUND_FLEXFN_H

#include "arcbounds.h"
#include <vector>

//Baseclass for FlexFN
class FlexFN : public ArcRemover {
 public:
  //Calculate field emission current I_r = alpha(r)*j(r)
  // [(#Superparticles/injection step) through annular area],
  // and add it to the array current
  void calcFN_current(double const Ez[], double* alpha, double* beta,
		      double* currentFN);
  //Inject superparticles into each ring, using the currentFN[]
  // (or a sum of those) calculated by calcFN_current()
  void injectFN(Particle pa[], size_t &np, double* currentFN, double const Ez[]);

};

//Ring emitter
class FlexFN_ring : public FlexFN {
 public:
  FlexFN_ring(std::vector<char*>& option);
  virtual ~FlexFN_ring();

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

  virtual const char* getName() const { return "FlexFN_ring"; };

 private:
  void initAlphaBeta();
  double* FN_alpha;
  double* FN_beta;
  double* FN_current;
  
  double alpha_ring, beta_ring;
  unsigned int idx1_ring, idx2_ring;
};

//Emitter summing over multiple distributions
class FlexFN_twoComp : public FlexFN {
 public:
  FlexFN_twoComp(std::vector<char*>& option);
  virtual ~FlexFN_twoComp();
  
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

  virtual const char* getName() const { return "FlexFN_twoComp"; };

 private:
  double** FN_alpha;
  double** FN_beta;
  double* FN_current;

  double alpha1, alpha2;
  double beta1, beta2;
  unsigned int idx1;
};


#endif
