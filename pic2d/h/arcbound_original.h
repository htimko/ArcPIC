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

  arcbounds_original.h:
  Particle boundary conditions,
  same as in Helga's thesis

***********************************************************************/

#ifndef ARCBOUND_ORIGINAL_H
#define ARCBOUND_ORIGINAL_H

#include "arcbounds.h"

#include <vector>

class ArcOriginal : public ArcRemover {
  //Implements the original arc boundary conditions from Helga's work
 public:
  ArcOriginal(std::vector<char*>& options);
  virtual ~ArcOriginal();

  void init(unsigned int nr, double zmin, double zmax, double rmax);
  void re_init(unsigned int nr, double zmin, double zmax, double rmax);

  virtual void print_par() const;

  virtual void inject_e(Particle pa[], size_t &np, double const Ez[]);
  virtual void inject_i(Particle pa[], size_t &np, double const Ez[], unsigned int sort) {}; //NOP
  virtual void inject_n(Particle pa[], size_t &np, double const Ez[]);

  virtual void remove_i(Particle pa[], size_t &np, unsigned int sort);
  virtual void remove_n(Particle pa[], size_t &np);
  //e- removal taken care of by ArcRemover

  virtual void timestep(unsigned int nstep, bool isOutputTimestep);

  //Write to the arcbounds.dat and arcbounds_original.dat output files
  virtual void writeFile_extras(unsigned int nstep) {
    writeFile_arcboundsOriginalDat(nstep);
  }
  
  //Save and restore backup data
  virtual void backup(FILE* file) {};
  virtual void restoreBackup(FILE* file) {};  

  virtual const char* getName() const { return "ArcOriginal"; };

 private:
  //Output file "arcbounds_original.dat"
  void writeFile_arcboundsOriginalDat(unsigned int nstep);
  void initFile_arcboundsOriginalDat();
  FILE* ofile_arcboundsOriginalDat;
  double E_tip_loc;

  //Settings from input file
  double beta_tip;      //Beta of the tip (erroded/melted over time)
  double beta_f;        //Flat surface
  double j_melt;        //Current where beta_tip becomes equal beta_flat

  bool   tipPi;         //Include the missing PI in the tip current
  double alpha_flat;    //Fractional of flat surface covered by emitters
  bool   fracInjectStep;//Inject electrons at fractional timestep (euler method)

  double Remission, Remission_theor; //Electron emission area
  
  double SEY;    //Single electron yield from copper
  
  double r_Cu_e; // ratio of copper evaporated to electrons FN-emitted
  
  double heatspike_threshold; //Threshold (in A/cm^2) for ion current through a cell
                              // on the cathode to trigger a heatspike
  double heatspike_yield;     //Yield in Cu neutrals / standard timestep
                              // (a standard timestep is approximately 1.77 fs,
                              // and arises at Omega_pe = 0.2, n_e = 4e18 cm^-3)
  
  //Needed data for sputtering
  std::vector<Sput> sput_cathode;
  std::vector<Sput> sput_cathode_SEY;
  std::vector<Sput> sput_anode;
  
  //Status for erosion & melting of tip
  bool has_melted;
  int erosion;
  int emitted;
  int emitted_tip;

  // 1..nr or 1e6*nr (no heatspike)
  double sigma_heatspike;

  //Injection velocities (based on thermal velocity)
  double v_inj_e;
  double v_inj_i;
  
  //Calculate sputtering for a single particle
  // Returns a sputtering object, where Y=r=0 if there is no sputtering
  inline Sput calc_sput(const Particle& pa, const double cs, double* current_enhancedY);
  //Inject netrals from sputtering on a single wall
  inline void inject_sput(Particle pa[], size_t &np,
			  std::vector<Sput> &sput, bool isCathode);
};


#endif
