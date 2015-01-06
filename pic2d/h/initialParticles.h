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

  initialParticles.h:
  Pre-fill the system with some particles
   (neutral background gas etc.)

  To implement your own initial particle distribution:
  1) Add the class declaration, inheriting InitialParticles
  2) Add it to the loader in InitialParticles::LoadInitialParticles()

***********************************************************************/

#ifndef INITIALPARTICLES_H
#define INITIALPARTICLES_H

#include "pic.h"

#include <stdio.h>
#include <vector>


class InitialParticles {
 public:
  //InitialParticles() {};
  //virtual ~InitialParticles() {};

  static InitialParticles* LoadInitialParticles(FILE* in_file);

  virtual void init() = 0;

  virtual void inject_e(Particle pa[], size_t &np) = 0;
  virtual void inject_n(Particle pa[], size_t &np) = 0;
  virtual void inject_i(Particle pa[], size_t &np, unsigned int sort) = 0;
  
  virtual void print_par() = 0;

  virtual const char* getName() const = 0;
};

class UniformRestricted : public InitialParticles {
 public:
  UniformRestricted (std::vector<char*> options);

  virtual void init();  

  virtual void inject_e(Particle pa[], size_t &np);
  virtual void inject_n(Particle pa[], size_t &np);
  virtual void inject_i(Particle pa[], size_t &np, unsigned int sort);
  
  virtual void print_par();

  virtual const char* getName() const {
    return "UniformRestricted";
  };

 private:
  double density; //Density in particles/cm^3
  double maxR; //Maximum R & Z coordinates in dz
  double maxZ;
  double minR; //Minimum ------ * * ------
  double minZ;
  
  bool doInject_e;
  bool doInject_i;
  bool doInject_n;

  double Tinj; //Injection temperature [eV]

  void injector(Particle pa[], size_t &np, double vInj);

  double vol;
  double Ldb;
  double N_sp;
  size_t num_inject;
};

#endif
