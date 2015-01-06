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

  arcbounds_original_newheatspike.h:
  Particle boundary conditions, same as in Helga's thesis,
  except for the heat spike sputtering part which has been updated,
  and tipPi which is locked to "yes"

***********************************************************************/

#ifndef ARCBOUND_ORIGINAL_NEWHEATSPIKE_H
#define ARCBOUND_ORIGINAL_NEWHEATSPIKE_H

#include "arcbounds.h"

#include <vector>

class ArcOriginalNewHS : public ArcRemover {
  //Implements the original arc boundary conditions from Helga's work
 public:
  ArcOriginalNewHS(std::vector<char*>& options);
  virtual ~ArcOriginalNewHS();

  void init(unsigned int nr, double zmin, double zmax, double rmax);
  void re_init(unsigned int nr, double zmin, double zmax, double rmax);

  virtual void print_par() const;

  virtual void inject_e(Particle pa[], size_t &np, double const Ez[]);
  virtual void inject_i(Particle pa[], size_t &np, double const Ez[], unsigned int sort) {}; //NOP
  virtual void inject_n(Particle pa[], size_t &np, double const Ez[]);

  virtual void remove_i(Particle pa[], size_t &np, unsigned int sort);
  virtual void remove_n(Particle pa[], size_t &np);
  //remove_e() taken care of by ArcRemover

  virtual void timestep(unsigned int nstep, bool isOutputTimestep);

  //Write to the arcbounds_original.dat output file
  virtual void writeFile_extras(unsigned int nstep) {
    writeFile_arcboundsOriginalDat(nstep);
  }
  
  //Save and restore backup data
  virtual void backup(FILE* file) {};
  virtual void restoreBackup(FILE* file) {};  

  virtual const char* getName() const { return "ArcOriginalNewHS"; };

 private:
  //Output file "arcbounds_original.dat"
  void writeFile_arcboundsOriginalDat(unsigned int nstep);
  FILE* ofile_arcboundsOriginalDat;
  double E_tip_loc;

  //Settings from input file
  double beta_tip;      //Beta of the tip (erroded/melted over time)
  double beta_f;        //Flat surface
  double j_melt;        //Current where beta_tip becomes equal beta_flat

  double alpha_flat;    //Fractional of flat surface covered by emitters
  bool   fracInjectStep;//Inject electrons at fractional timestep (euler method)

  double Remission, Remission_theor; //Electron emission area  
  double Rborder;                    //Inner radius of flat emission area

  double SEY;    //Single electron yield from copper
  
  double r_Cu_e;         // Ratio of copper evaporated to electrons FN-emitted
  double r_Cu_e_flat;    // Ratio of copper evaporated to electrons FN-emitted from flat surface
  bool evap_flat_center; // Evaporate copper from flat surface inside Remission?

  double heatspike_threshold; //Threshold (in particles/cm^2/s) for Cu and Cu+ flux through a cell
                              // on the cathode to trigger a heatspike
  double heatspike_yield_ratio; //Yield ratio "p" -- number of neutrals emitted / number of impactors

  char heatspike_model; //a: Activate if one or more cells are above threshold.
                        //   Count impactors in these cells only.
                        //b: Active as for a, but count impactors also in cell inside activated cells.
                        //c: Activate if average current though some circular area bounded by a cell outer edge
                        //   is above threshold. Count impactors in this area.

  //Needed data for sputtering
  std::vector<Sput> sput_cathode;
  std::vector<Sput> sput_cathode_SEY;
  std::vector<Sput> sput_anode;
  unsigned int* sput_cathode_current; //For determining heatspike

  //Status for erosion & melting of tip
  int erosion;
  bool do_erosion;
  int emitted_tip_melting;
  bool has_melted;
  //For evaporation [number of electron superparticles]
  int emitted_tip_evap;
  int* emitted_flat_evap;
  //For output [number of electron superparticles]
  int emitted_tip_output;
  int emitted_flat_output;
  int emitted_SEY_output;
  //[number of neutral superparticles]
  int emitted_evap_output;
  int emitted_sputter_cat_output; //cathode sputtering
  int emitted_sputter_ano_output; //  anode sputtering
  int emitted_heatspike_output;   //heatspike sputtering (on cathode)
  
  // 1..nr or 0 (no heatspike)
  unsigned int heatspike_sigma;
  unsigned int heatspike_incident; // [number of neutral superparticles incident]

  //Injection velocities (based on thermal velocity)
  double v_inj_e;
  double v_inj_i;

  //Calculate sputtering for a single particle
  // Returns a sputtering object, where Y=r=0 if there is no sputtering
  inline Sput calc_sput(const Particle& pa, const double cs, unsigned int* current_enhancedY);
  //Inject netrals from sputtering on a single wall.
  // The sputtering array is automatically reset.
  inline void inject_sput(Particle pa[], size_t &np,
			  std::vector<Sput> &sput, bool isCathode);  
};


#endif
