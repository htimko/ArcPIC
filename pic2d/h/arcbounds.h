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

  arcbounds.h:
  Particle boundary condition classes

***********************************************************************/

#ifndef ARCBOUNDS_H
#define ARCBOUNDS_H

#include "pic.h"
#include "dim.h"

#include <vector>
#include <iostream>

#include <cstdlib>
#include <cstdio>

class ArcBounds {
 public:
  //You should also provide your own constructor to deal with input data
  // (this constructor will anyway be called)
  ArcBounds();
  virtual ~ArcBounds();

  //This function identifes and loads the right ArcBounds class
  static ArcBounds* LoadArcBounds(FILE* in_file);
  
  //Called by calc_parameters_2D(), setups rescaling etc.
  virtual void init(unsigned int nr, double zmin, double zmax, double rmax);
  //Called when re-initializing from backup
  virtual void re_init(unsigned int nr, double zmin, double zmax, double rmax);
  //Called by print_par_2D()
  virtual void print_par() const = 0;
  
  //Save and restore backup data
  virtual void backup(FILE* file)        = 0;
  virtual void restoreBackup(FILE* file) = 0;
  
  //Remove out-of-bounds particles
  virtual void remove_e(Particle pa[], size_t &np) = 0;
  virtual void remove_i(Particle pa[], size_t &np, unsigned int sort) = 0;
  virtual void remove_n(Particle pa[], size_t &np) = 0;
  
  //Inject particles
  virtual void inject_e(Particle pa[], size_t &np, double const Ez[]) = 0;
  virtual void inject_i(Particle pa[], size_t &np, double const Ez[], unsigned int sort) = 0;
  virtual void inject_n(Particle pa[], size_t &np, double const Ez[]) = 0;
  
  //Get transported charge, advance to next timestep
  virtual const double getDeltaQ() const;
  virtual void timestep(unsigned int nstep, bool isOutputTimestep);
  
  virtual const char* getName() const = 0;

  inline const unsigned int getFileTimestepHist() const {return file_timestep_hist;};

 protected:
  //Hook for inheriting classes to make their own extra output files
  virtual void writeFile_extras(unsigned int nstep) {};

  //Main output file "arcbounds.dat",
  // with lots of data about added/removed particles
  virtual void writeFile_arcboundsDat(unsigned int nstep);
  unsigned int file_timestep;
  //Arrays for keeping track of particle injection/removal
  // These values are always positive, but declared as signed
  // to avoid problems when subtracting them
  // [0] : z<Zmin "wall1" (cathode)
  // [1] : z>Zmax "wall2" (anode)
  // [2] : r>Rmax "Infinity";
  //How many particles where added at each wall
  // (nothing at infinity)
  int injected_e[2];
  int injected_i[NSpecies-1][2];
  int injected_n[2];
  //How many particles where removed at each wall  
  int removed_e[3];
  int removed_i[NSpecies-1][3];
  int removed_n[3];
  //Total current density of each species at each wall
  // (nothing at infinity)
  int current_e[2];
  int current_i[NSpecies-1][2];
  int current_n[2];
  //Reset to zero
  void resetCountingArrays();

  //Current density as function of r, binned same as mesh
  void writeFile_currentHistos(unsigned int nstep, bool isOutputTimestep);
  //CurrentHistos are written at timestep modulo this value.
  // If set to 0, only write on common output timestep.
  unsigned int file_timestep_hist;
  int* current_cathode;
  int* current_anode;

  //Dump all particles which are removed
  void writeFile_removedParticles(unsigned int nstep);
 protected:
  //Removed particle dumps for writeFile_removedParticles()
  std::vector<Particle> removedElectrons;
  std::vector<Particle> removedIons;
  std::vector<Particle> removedNeutrals;

  //System size
  unsigned int nr;
  double zmin, zmax, rmax;
 private:
  FILE* ofile_arcboundsDat;
  
  FILE* ofile_currentHist_cathode;
  FILE* ofile_currentHist_anode;
  
  FILE* ofile_removedParticles_electrons;
  FILE* ofile_removedParticles_ions;
  FILE* ofile_removedParticles_neutrals;
};

class ArcRemover : public ArcBounds {
  // Pure virtual class removing particles at the boundaries and counting them,
  // no sputtering etc.
 
  virtual void remove_e(Particle pa[], size_t &np);
  virtual void remove_i(Particle pa[], size_t &np, unsigned int sort);
  virtual void remove_n(Particle pa[], size_t &np);
  
 protected:
  //This function does the actual work
  void remover(Particle pa[], size_t &np,
	       int removed[], int current[],
	       int chargeSign,
	       std::vector<Particle>& removedVector );
};

class ArcDummy : public ArcBounds {
  //No particles injected or removed - but not pure virtual.
  // Usefull for testing or as a template
  
 public:
  ArcDummy(std::vector<char*>& options);

  virtual void print_par() const {};
  
  //Save and restore backup data
  virtual void backup(FILE* file) {};
  virtual void restoreBackup(FILE* file) {};

  //Remove out-of-bounds particles
  virtual void remove_e(Particle pa[], size_t &np);
  virtual void remove_i(Particle pa[], size_t &np, unsigned int sort);
  virtual void remove_n(Particle pa[], size_t &np);  

  //Inject particles
  virtual void inject_e(Particle pa[], size_t &np, double const Ez[]);
  virtual void inject_i(Particle pa[], size_t &np, double const Ez[], unsigned int sort);
  virtual void inject_n(Particle pa[], size_t &np, double const Ez[]);
  
  virtual const char* getName() const { return "ArcDummy"; };
};

class ArcSimple : public ArcRemover {
 public:
  //Very simple ArcBound; only inject electrons at a steady, setable rate.
  ArcSimple(std::vector<char*>& options);
  virtual void print_par() const;
  
  //Inject particles
  virtual void inject_e(Particle pa[], size_t& np, double const Ez[]);
  virtual void inject_i(Particle pa[], size_t &np, double const Ez[], unsigned int sort);
  virtual void inject_n(Particle pa[], size_t &np, double const Ez[]);

  
  //Save and restore backup data
  virtual void backup(FILE* file) {};
  virtual void restoreBackup(FILE* file) {};  

  virtual const char* getName() const { return "ArcSimple"; };

 private:
  unsigned int ne_inject, ni_inject, nn_inject;

  void injector(Particle pa[], size_t &np, unsigned int n_inject);
};



#endif
