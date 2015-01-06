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

  arcbounds_flexFN_origNeutrals.h:
  Particle boundary conditions,
  Fowler-Nordheim with beta and alpha functions of r,
  with sputtering and evaporation similar to Helga's thesis.

***********************************************************************/

#ifndef ARCBOUND_FLEXFN_ORIGNEUTRALS_H
#define ARCBOUND_FLEXFN_ORIGNEUTRALS_H

#include "arcbounds.h"
#include "arcbound_flexFN.h"
#include <vector>

//FlexFN-type emitter summing over multiple distributions,
// using original neutral sputtering scheme
// and per-cell neutral evaproation
class FlexFN_twoComp_origNeutrals : public FlexFN {
 public:
  FlexFN_twoComp_origNeutrals(std::vector<char*>& option);
  virtual ~FlexFN_twoComp_origNeutrals();
  
  virtual void init   (unsigned int nr, double zmin, double zmax, double rmax);
  virtual void re_init(unsigned int nr, double zmin, double zmax, double rmax);

  virtual void print_par() const;

  virtual void inject_e(Particle pa[], size_t &np, double const Ez[]);
  virtual void inject_i(Particle pa[], size_t &np, double const Ez[],
			unsigned int sort) {}; //NOP
  virtual void inject_n(Particle pa[], size_t &np, double const Ez[]);

  virtual void remove_i(Particle pa[], size_t &np, unsigned int sort);
  virtual void remove_n(Particle pa[], size_t &np);

  //Save and restore backup data
  virtual void backup(FILE* file) {};
  virtual void restoreBackup(FILE* file) {};  

  virtual const char* getName() const { return "FlexFN_twoComp_origNeutrals"; };

 private:
  //Fowler-Nordheim
  void initAlphaBeta();

  double** FN_alpha;
  double** FN_beta;
  double* FN_current;

  double alpha1, alpha2;
  double beta1, beta2;
  unsigned int idx1;
  
  //SEY on cathode
  double SEY;    //Single electron yield from copper
  std::vector<Sput> sput_cathode_SEY;

  //Neutral injection/Evaporation
  double r_Cu_e; // ratio of copper evaporated to electrons FN-emitted
  double* FN_current_sum; //Sum FN_current over neutral injection timestep
                          // in order to find how many neutrals should be injected.

  //Needed data for sputtering
  std::vector<Sput> sput_cathode;
  std::vector<Sput> sput_anode;
  
  bool doHeatspike;

  //Calculate sputtering for a single particle
  // Returns a sputtering object, where Y=r=0 if there is no sputtering
  inline Sput calc_sput(const Particle& pa, const double cs, double* current_enhancedY);
  //Inject netrals from sputtering on a single wall
  inline void inject_sput(Particle pa[], size_t &np,
			  std::vector<Sput> &sput,
			  bool isCathode, double v_inj);

};

#endif
