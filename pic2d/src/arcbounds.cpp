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
  
  arcbounds.cpp:
  Defines the common services and default function implementations
  for ArcBounds classes. These classes implement surface physics models,
  and which class to load can be selected by the input file.
  Some classes inheriting from and implementing all or part of the virtual methods
  in ArcBounds are also defined in this file, such as ArcRemover, ArcDummy and
  ArcSimple.

***********************************************************************/

#include <iostream>
#include <vector>

#include <cstdlib>
#include <cstring>
#include <cmath>

#include "arcbounds.h"
#include "arcbound_original.h"
#include "arcbound_original_newheatspike.h"
#include "arcbound_flexFN.h"
#include "arcbound_flexFN_origNeutrals.h"
#include "arcbound_flexFN2.h"
#include "arcbound_flexFN2_tableNeutrals.h"
#include "arcbound_sumTipFNField.h"

#include "input.h"

#include "random.h"

#define   XTRN extern
#include "var.h"
#include "outp.h" //Filenames
#undef XTRN

using namespace std;

// ******** Implementation of ArcBounds ******************

//Static functions in class ArcBounds
ArcBounds* ArcBounds::LoadArcBounds(FILE* in_file) {
  vector<char*> options;
  char* arcBoundsType = readInputSection(in_file, options);

  ArcBounds* ret = NULL;
  
  if (! strncmp(arcBoundsType, "ArcDummy", NAME_MAXLEN) ) {
    ret = new ArcDummy(options);
  }
  else if (! strncmp(arcBoundsType, "ArcSimple", NAME_MAXLEN) ) {
    ret = new ArcSimple(options);
  }
  else if (! strncmp(arcBoundsType, "ArcOriginal", NAME_MAXLEN) ) {
    ret = new ArcOriginal(options);
  }
  else if (! strncmp(arcBoundsType, "ArcOriginalNewHS", NAME_MAXLEN) ) {
    ret = new ArcOriginalNewHS(options);
  }
  else if (! strncmp(arcBoundsType, "FlexFN_ring", NAME_MAXLEN) ) {
    ret = new FlexFN_ring(options);
  }
  else if (! strncmp(arcBoundsType, "FlexFN_twoComp", NAME_MAXLEN) ) {
    ret = new FlexFN_twoComp(options);
  }
  else if (! strncmp(arcBoundsType, "FlexFN_twoComp_origNeutrals", NAME_MAXLEN) ) {
    ret = new FlexFN_twoComp_origNeutrals(options);
  }
  else if (! strncmp(arcBoundsType, "FlexFN2_ring", NAME_MAXLEN) ) {
    ret = new FlexFN2_ring(options);
  }
  else if (! strncmp(arcBoundsType, "FlexFN2_tableNeutrals", NAME_MAXLEN) ) {
    ret = new FlexFN2_tableNeutrals(options);
  }
  else if (! strncmp(arcBoundsType, "SumTipFNField", NAME_MAXLEN) ) {
    ret = new SumTipFNField(options);
  }

  else {
    cout << "Unknown ArcBounds type '" << arcBoundsType << "'" << endl;
    exit(1);
  }

  //Cleanup
  delete[] arcBoundsType;
  for (size_t i = 0; i < options.size(); i++) {
    delete[] options[i];
  }
  options.clear();
  
  return ret;
}

ArcBounds::ArcBounds() : file_timestep(1), file_timestep_hist(0), ofile_arcboundsDat(NULL), ofile_currentHist_cathode(NULL), ofile_currentHist_anode(NULL), ofile_removedParticles_electrons(NULL), ofile_removedParticles_ions(NULL), ofile_removedParticles_neutrals(NULL) {
  resetCountingArrays();
}
ArcBounds::~ArcBounds() {
  if (ofile_arcboundsDat != NULL) fclose(ofile_arcboundsDat);

  delete[] current_cathode;
  delete[] current_anode;
  fclose(ofile_currentHist_cathode);
  fclose(ofile_currentHist_anode);
}

void ArcBounds::init(unsigned int nr, double zmin, double zmax, double rmax) {
  this->nr   = nr;
  this->zmin = zmin;
  this->zmax = zmax;
  this->rmax = rmax;

  this->current_cathode = new int[this->nr];
  this->current_anode = new int[this->nr];
  for (unsigned int i = 0; i < nr; i++){
    current_cathode[i] = current_anode[i] = 0;
  }
}
void ArcBounds::re_init(unsigned int nr, double zmin, double zmax, double rmax) {
  init(nr, zmin, zmax, rmax);
}

void ArcBounds::writeFile_arcboundsDat(unsigned int nstep) {
  if (ofile_arcboundsDat == NULL) {
    ofile_arcboundsDat = fopen("arcbounds.dat", "w");
    fprintf(ofile_arcboundsDat, "## ArcBounds : %s ##\n", this->getName());
    fprintf(ofile_arcboundsDat, "## If double index: First is species (H, Cu), 2nd is wall ID (cathode, anode, rmax). \n");
    fprintf(ofile_arcboundsDat, "## If single index: Only wall ID. \n");
    fprintf(ofile_arcboundsDat, "## timestep injected_e[0,1] injected_i[0,...,%u][0,1] injected_n[0,1] removed_e[0,1,2] removed_i[0,...%u][0,1,2] removed_n[0,1,2] current_e[0,1] current_i[0,..%u][0,1] current_n[0,1]\n", NSpecies-2, NSpecies-2, NSpecies-2);
    fflush(ofile_arcboundsDat);
  }
  
  if ( nstep % file_timestep == 0) {
    fprintf(ofile_arcboundsDat, "%10d %6d %6d ", nstep, injected_e[0], injected_e[1]);
    for (size_t i = 0; i < NSpecies-1; i++) {
      fprintf(ofile_arcboundsDat, "%6d %6d ", injected_i[i][0], injected_i[i][1]);
    }
    fprintf(ofile_arcboundsDat, "%6d %6d %6d %6d %6d ", injected_n[0], injected_n[1], removed_e[0], removed_e[1], removed_e[2]);
    for (size_t i = 0; i < NSpecies-1; i++) {
      fprintf(ofile_arcboundsDat, "%6d %6d %6d ", removed_i[i][0], removed_i[i][1], removed_i[i][2]);
    }
    fprintf(ofile_arcboundsDat, "%6d %6d %6d %6d %6d ", removed_n[0], removed_n[1], removed_n[2],
	    current_e[0], current_e[1]);
    for (size_t i = 0; i < NSpecies-1; i++) {
      fprintf(ofile_arcboundsDat, "%6d %6d ", current_i[i][0], current_i[i][1]);
    }
    fprintf(ofile_arcboundsDat, "%6d %6d \n ", current_n[0], current_n[1]);
    
    fflush(ofile_arcboundsDat);
    
    resetCountingArrays();
  }
}
void ArcBounds::resetCountingArrays() {
  injected_e[0]=injected_e[1] = 0;
  injected_n[0]=injected_n[1] = 0;
  removed_e[0]=removed_e[1]=removed_e[2] = 0;
  removed_n[0]=removed_n[1]=removed_n[2] = 0;
  current_e[0]=current_e[1] = 0;
  current_n[0]=current_n[1] = 0;
  for (size_t i = 0; i<NSpecies-1; i++) {
    injected_i[i][0]=injected_i[i][1] = 0;
    removed_i[i][0]=removed_i[i][1]=removed_i[i][2] = 0;
    current_i[i][0]=current_i[i][1] = 0;
  }
}

void ArcBounds::writeFile_currentHistos(unsigned int nstep, bool isOutputTimestep) {
  if (! isOutputTimestep ) {
    if (file_timestep_hist == 0) return; //In this case only activate on engage
    if (nstep % file_timestep_hist != 0) return;
  }

  //Write for the cathode  
  if (ofile_currentHist_cathode == NULL) {
    ofile_currentHist_cathode = fopen("out/jhist_cathode.dat", "w");
    
    fprintf(ofile_currentHist_cathode, "## ArcBounds : %s ##\n", this->getName());
    fprintf(ofile_currentHist_cathode, "## Cathode, nr=%u, rmax=%g ## \n", this->nr, this->rmax);
    fprintf(ofile_currentHist_cathode, "## timestep I[0] I[1] I[2] ... I[nr-1] [num. superparticles / output step, pos. dir = 'electrons anode->cathode']\n");
    fflush(ofile_currentHist_cathode);
  }
  fprintf(ofile_currentHist_cathode, "%10d ", nstep);
  for (unsigned int i = 0; i < nr; i++) {
    fprintf(ofile_currentHist_cathode, "%d ", current_cathode[i]);
    current_cathode[i] = 0;
  }
  fprintf(ofile_currentHist_cathode, "\n");
  fflush(ofile_currentHist_cathode);

  //Write for the anode
  if (ofile_currentHist_anode == NULL) {
    ofile_currentHist_anode = fopen("out/jhist_anode.dat", "w");
    
    fprintf(ofile_currentHist_anode, "## ArcBounds : %s ##\n", this->getName());
    fprintf(ofile_currentHist_anode, "## Anode, nr=%u, rmax=%g ## \n", this->nr, this->rmax);
    fprintf(ofile_currentHist_anode, "## timestep I[0] I[1] I[2] ... I[nr-1] [num. superparticles / output step, pos. dir = 'electrons anode->cathode']\n");
    fflush(ofile_currentHist_anode);
  }
  fprintf(ofile_currentHist_anode, "%10d ", nstep);
  for (unsigned int i = 0; i < nr; i++) {
    fprintf(ofile_currentHist_anode, "%d ", current_anode[i]);
    current_anode[i] = 0;
  }
  fprintf(ofile_currentHist_anode, "\n");
  fflush(ofile_currentHist_anode);
  
}

void ArcBounds::writeFile_removedParticles(unsigned int nstep) {
  if (ofile_removedParticles_electrons == NULL) {
    ofile_removedParticles_electrons = fopen("out/removed_e.dat", "w");
    fprintf(ofile_removedParticles_electrons, "## timestep z r vz vr vt\n");
    ofile_removedParticles_ions      = fopen("out/removed_i.dat", "w");
    fprintf(ofile_removedParticles_ions, "## timestep z r vz vr vt\n");
    ofile_removedParticles_neutrals  = fopen("out/removed_n.dat", "w");
    fprintf(ofile_removedParticles_neutrals, "## timestep z r vz vr vt\n");
  }
  
  for (size_t i = 0; i < removedElectrons.size(); i++) {
    PhasespaceVar& p = removedElectrons[i].p;
    fprintf(ofile_removedParticles_electrons, "%10d %.5e %.5e %.5e %.5e %.5e\n", nstep, p.z, p.r, p.vz, p.vr, p.vt);
  }
  fflush(ofile_removedParticles_electrons);
  removedElectrons.clear();

  for (size_t i = 0; i < removedIons.size(); i++) {
    PhasespaceVar& p = removedIons[i].p;
    fprintf(ofile_removedParticles_ions, "%10d %.5e %.5e %.5e %.5e %.5e\n", nstep, p.z, p.r, p.vz, p.vr, p.vt);
  }
  fflush(ofile_removedParticles_ions);
  removedIons.clear();

  for (size_t i = 0; i < removedNeutrals.size(); i++) {
    PhasespaceVar& p = removedNeutrals[i].p;
    fprintf(ofile_removedParticles_neutrals, "%10d %.5e %.5e %.5e %.5e %.5e\n", nstep, p.z, p.r, p.vz, p.vr, p.vt);
  }
  fflush(ofile_removedParticles_neutrals);
  removedNeutrals.clear();
}

const double ArcBounds::getDeltaQ() const {
  //Simple implementation:
  // Transported charge this timestep = charge transport through cathode
  // Would be better to use Shockley-Ramo.
  double deltaQ = injected_e[0]-removed_e[0];
  for (size_t i = 0; i<NSpecies-1; i++) {
    deltaQ += removed_i[i][0]-injected_i[i][0];
  }
  return deltaQ;
}

void ArcBounds::timestep(unsigned int nstep, bool isOutputTimestep) {
  //Kick off file writing
  writeFile_arcboundsDat(nstep);
  writeFile_removedParticles(nstep);
  writeFile_currentHistos(nstep, isOutputTimestep);
  writeFile_extras(nstep);
}

// ******** Implementation of ArcRemover ******************
void ArcRemover::remove_e(Particle pa[], size_t &np) {
  remover(pa, np, removed_e, current_e, -1, removedElectrons);
}

void ArcRemover::remove_i(Particle pa[], size_t &np, unsigned int sort) {
  remover(pa, np, removed_i[sort], current_i[sort], 1, removedIons);
}

void ArcRemover::remove_n(Particle pa[], size_t &np) {
  remover(pa, np, removed_n, current_n, 0, removedNeutrals);
}

void ArcRemover::remover(Particle pa[], size_t &np,
			 int removed[],
			 int current[],
			 int chargeSign,
			 std::vector<Particle>& removedVector) {
  int jr;
  size_t n_lost = 0;
  
  for (size_t n=0; n<np; n++ ) {
    if ( pa[n].p.r >= rmax ) {
      removed[2]++;
      removedVector.push_back(pa[n]);
      n_lost++;
      continue; 
    }
    else if ( pa[n].p.z >= zmax ) {
      removed[1]++;     
      jr = int (pa[n].p.r);
      current[1] += 2*jr+1;
      current_anode[jr] -= chargeSign;
      removedVector.push_back(pa[n]);
      n_lost++;
      continue; 
    }
    else if ( pa[n].p.z < zmin ) {
      removed[0]++;
      jr = int (pa[n].p.r);
      current[0] += 2*jr+1;
      current_cathode[jr] += chargeSign;
      removedVector.push_back(pa[n]);
      n_lost++;
      continue; 
    }
    //Implicit else: keep this particle
    pa[n-n_lost].p = pa[n].p;
  }
  
  np -= n_lost;
}

// ******** Implementation of ArcDummy ******************
ArcDummy::ArcDummy(std::vector<char*>& options) {
  cout << "In ArcDummy constructor" << endl;
  if (options.size() != 0) {
    cout << "Error in ArcDummy(): Expected 0 options, got "
	 << options.size() << endl;
    exit(1);
  }
}

void ArcDummy::remove_e(Particle pa[], size_t &np) {
  cout << "remove_e" << endl;
}
void ArcDummy::remove_i(Particle pa[], size_t &np, unsigned int sort) {
  cout << "remove_i, sort=" << sort << endl;
}
void ArcDummy::remove_n(Particle pa[], size_t &np) {
  cout << "remove_n" << endl;
}

void ArcDummy::inject_e(Particle pa[], size_t &np, double const Ez[]) {
  cout << "inject_e" << endl;
}
void ArcDummy::inject_i(Particle pa[], size_t &np, double const Ez[], unsigned int sort) {
  cout << "inject_i, sort=" << sort << endl;
}
void ArcDummy::inject_n(Particle pa[], size_t &np, double const Ez[]) {
  cout << "inject_n" << endl;
}

// ******** Implementation of ArcSimple ******************
ArcSimple::ArcSimple(vector<char*>& options) {
  if (options.size() != 5) {
    cout << "Error in ArcSimple(): Expected 5 options, got "
	 << options.size() << endl;
    exit(1);
  }

  sscanf(options[0], "%*[^:]%*[:] %u", &(this->ne_inject));
  sscanf(options[1], "%*[^:]%*[:] %u", &(this->ni_inject));
  sscanf(options[2], "%*[^:]%*[:] %u", &(this->nn_inject));
  sscanf(options[3], "%*[^:]%*[:] %u", &(this->file_timestep));
  sscanf(options[4], "%*[^:]%*[:] %u", &(this->file_timestep_hist));
}
void ArcSimple::print_par() const {
  printf( " - ne_inject            %u \n", this->ne_inject );
  printf( " - ni_inject            %u \n", this->ni_inject );
  printf( " - nn_inject            %u \n", this->nn_inject );
  printf( " - file_timestep        %u \n", this->file_timestep );
  printf( " - file_timestep_hist   %u \n", this->file_timestep_hist );
}

void ArcSimple::inject_e(Particle pa[], size_t &np, double const Ez[]) {
  injector(pa,np,ne_inject);
  injected_e[0] += ne_inject;
  current_e [0] += ne_inject; //All injected in cell (0,0)
  current_cathode[0] += ne_inject;
}
void ArcSimple::inject_i(Particle pa[], size_t &np, double const Ez[], unsigned int sort) {
  if ( sort != 1) return; //Only make Cu
  injector(pa,np,ni_inject);
  injected_i[1][0] += ni_inject;
  current_i [1][0] += ni_inject; //All injected in cell (0,0)
  current_cathode[0] -= ni_inject;
}
void ArcSimple::inject_n(Particle pa[], size_t &np, double const Ez[]) {
  injector(pa,np,nn_inject);
  injected_n[0] += nn_inject;
}

void ArcSimple::injector(Particle pa[], size_t &np, unsigned int n_inject) {
  //Insert particles with zero velocity randomly into first cell
  if (np + n_inject >= NPART ) {
    cout << "Error in ArcSimple::injector(): Particle overflow" << endl;
    exit(1);
  }

  for (unsigned int i = 0; i < n_inject; i++) {
    pa[np+i].p.r = RAND*dr;
    pa[np+i].p.z = RAND*dz;
    
    pa[np+i].p.vr = pa[np+i].p.vt = pa[np+i].p.vz = 0.0;
    pa[np+i].p.m = 1;
  }
  
  np += n_inject;
}

