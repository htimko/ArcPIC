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

  arcbounds_flexFN2_tableNeutrals.cpp:
  Particle boundary conditions,
  Neutral sputtering from tabulated distrubutions,
  and FlexFN2 electron injection.

***********************************************************************/

#include "arcbound_flexFN2_tableNeutrals.h"
#include "mydef.h"
#include "random.h"

#define   XTRN extern
#include "var.h"
#include "arrays1.h"
#undef XTRN

#include <cmath>

using namespace std;


// ******** Implementation of FlexFN2_tableNeutrals ******************
FlexFN2_tableNeutrals::FlexFN2_tableNeutrals(std::vector<char*>& options) {
  if (options.size() != 6) {
    cout << "Error in FlexFN2_tableNeutrals(): Expected 6 options, got "
	 << options.size() << endl;
    exit(1);
  }

  sscanf(options[0], "%*[^:]%*[:] %lg", &(this->emitter_radius));
  sscanf(options[1], "%*[^:]%*[:] %lg", &(this->emitter_alpha));
  sscanf(options[2], "%*[^:]%*[:] %lg", &(this->emitter_beta));

  sscanf(options[3], "%*[^:]%*[:] %lg", &(this->flat_alpha));
  sscanf(options[4], "%*[^:]%*[:] %lg", &(this->flat_beta));
  
  sscanf(options[5], "%*[^:]%*[:] %u",  &(this->file_timestep));
}
void FlexFN2_tableNeutrals::print_par() const {
  printf( " - emitter_radius                %g \n", emitter_radius);
  printf( " - emitter_alpha                 %g \n", emitter_alpha);
  printf( " - emitter_beta                  %g \n", emitter_beta);
  printf( " - flat_alpha                    %g \n", flat_alpha);
  printf( " - flat_beta                     %g \n", flat_beta);
  printf( " - file_timestep                 %u \n", file_timestep );
}

void FlexFN2_tableNeutrals::init(unsigned int nr, double zmin, double zmax, double rmax) {
  ArcBounds::init(nr,zmin,zmax,rmax);
}
void FlexFN2_tableNeutrals::re_init(unsigned int nr, double zmin, double zmax, double rmax) {
  FlexFN2_tableNeutrals::init(nr, zmin, zmax, rmax);
}

void FlexFN2_tableNeutrals::inject_e(Particle pa[], size_t &np, double const Ez[]) {
  injectFNring(pa, np, emitter_alpha, emitter_beta, Ez, 0, emitter_radius);
  injectFNring(pa, np, flat_alpha,    flat_beta,    Ez, 0, nr);
}

void FlexFN2_tableNeutrals::inject_n(Particle pa[], size_t &np, double const Ez[]) {
  //Here goes the sputtering injection. Probably should also overwrite remove_i and remove_n
}
