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

  circuit.h:
  Models for the external circuit.

  To implement your own circuit model class:
  1) Add the class declaration, inheriting Circuit
  2) Add it to the loader in Circuit::LoadCircuit()

***********************************************************************/

#ifndef CIRCUIT_H
#define CIRCUIT_H

#include <cstdio>
#include <vector>

class Circuit {
  /*
   * Interface specification for an external circuit.
   * This class i virtual, meaning that you can't instanciate an object from it;
   * its just here to be inherited from.
   */
 public:
  Circuit() : ofile(NULL), file_timestep(1), deltaQ_accumulated(0.0) { };
  virtual ~Circuit() {fclose(ofile);};
  
  //Called by calc_parameters_2D(), does rescaling etc.
  virtual void init()      = 0;
  //Called when re-initializing from backup
  virtual void re_init()   = 0;
  //Called by print_par_2D()
  virtual void print_par() = 0;
  
  //Calculate the voltage on the gap in the next timestep
  virtual void timestep(double deltaQ, unsigned int nstep) = 0;
  //Write to the circuit.dat output file
  virtual void initFile();
  virtual void writeFile(unsigned int nstep);

  //Functions to actually get the new voltages
  inline const double getU0() const { return this->U0;  };
  inline const double getUNz() const { return this->UNz; };
  
  //Save and restore backup data
  virtual void backup(FILE* file)        = 0;
  virtual void restoreBackup(FILE* file) = 0;
  
  //Which circuit is this really?
  virtual const char* getName() const = 0;
  
  //This function identifes and loads the right circuit class
  static Circuit* LoadCircuit(FILE* in_file);
  
 protected: // These are accessed by child classes
  double U0;
  double UNz;
  
  FILE* ofile;
  unsigned int file_timestep;
  double deltaQ_accumulated;
};

class FixedVoltage : public Circuit {
  /*
   * Very simple external circuit: Just keep a fixed voltage on the gap.
   */
 public:
  FixedVoltage( std::vector<char*> options );
  
  virtual void init();
  virtual void re_init() { ofile = fopen("circuit.dat", "a"); init(); };
  virtual void print_par();
  virtual void timestep(double deltaQ, unsigned int nstep);
  
  //These are not needed in this function
  virtual void backup(FILE* file) {};
  virtual void restoreBackup(FILE* file) {};

  virtual const char* getName() const { return "FixedVoltage"; };

};

class FixedVoltage_resistor : public Circuit {
  /*
   * Fixed supply voltage and a series resistor
   */
 public:
  FixedVoltage_resistor( std::vector<char*> options );
  
  virtual void init();
  virtual void re_init() { ofile = fopen("circuit.dat", "a"); init(); };
  virtual void print_par();
  
  virtual void timestep(double deltaQ, unsigned int nstep);
  
  //These are not needed in this function
  virtual void backup(FILE* file) {};
  virtual void restoreBackup(FILE* file) {};

  virtual const char* getName() const { return "FixedVoltage_resistor"; };
  
 private:
  double Vin;     //Input voltage
  double Rseries, Rseries_dim; //Series resistance
};

class FixedVoltage_resistorCapacitor : public Circuit {
  /*
   * Fixed supply voltage, a series resistor, and gap capacitance.
   */
 public:
  FixedVoltage_resistorCapacitor( std::vector<char*> options );
  
  virtual void init();
  virtual void re_init() { ofile = fopen("circuit.dat", "a"); init(); };
  virtual void print_par();

  virtual void initFile();
  virtual void writeFile(unsigned int nstep);

  virtual void timestep(double deltaQ, unsigned int nstep);
  
  //These are not needed in this function
  virtual void backup(FILE* file) {};
  virtual void restoreBackup(FILE* file) {};

  virtual const char* getName() const { return "FixedVoltage_resistorCapacitor"; };
  
 private:
  double Vin;                  // Input voltage
  double Rseries, Rseries_dim; // Series resistance
  double Cgap, Cgap_dim;       // Gap capacitance
  double I_circ;               // Circuit current [superparticles/timestep]
};

class FixedVoltage_ramp : public Circuit {
  /*
   * Same as fixedVoltage, but ramps the voltage linearly
   * over rampTime [ns] until UNz-U0 is reached
   */
 public:
  FixedVoltage_ramp( std::vector<char*> options );
  
  virtual void init();
  virtual void re_init() { ofile = fopen("circuit.dat", "a"); init(); };
  virtual void print_par();
  
  virtual void timestep(double deltaQ, unsigned int nstep);
  
  //These are not needed in this function
  virtual void backup(FILE* file) {};
  virtual void restoreBackup(FILE* file) {};

  virtual const char* getName() const { return "FixedVoltage_ramp"; };
  
 private:
  double Vfinal;
  double rampTime;
  unsigned int rampTime_steps;
};

class TwoCapacitors : public Circuit {
  /*
   * This is the circuit from Helga's thesis
   * (with some modifications)
   */
 public:
  TwoCapacitors( std::vector<char*> options );

  virtual void init();
  virtual void re_init();
  virtual void print_par();

  virtual void initFile();
  virtual void writeFile(unsigned int nstep);

  virtual void timestep(double deltaQ_C, unsigned int nstep);
  
  virtual void backup(FILE* file);
  virtual void restoreBackup(FILE* file);

  virtual const char* getName() const { return "TwoCapacitors"; };
  
 private:
  double U0_dim;
  double UNz_dim;
  double C_ext , C_ext_dim;
  double C_ext2, C_ext2_dim;
  double R_ext , R_ext_dim; 
  double Q_ext;

  //TODO: Add these to backup
  double t_change_dim;
  unsigned int t_change;
  
  bool avoid_negative;
  
  bool check_C;
};

class FrequencyExcitation : public Circuit {
  /*
   * Excite anode with V(t) = UNz*sin(2*pi*f0 * t)
   * Cicuit current does NOT include charging of electrodes.
   */

 public:  
  FrequencyExcitation( std::vector<char*> options );
  
  virtual void init();
  virtual void re_init() { ofile = fopen("circuit.dat", "a"); init(); };
  virtual void print_par();
  virtual void timestep(double deltaQ, unsigned int nstep);
  
  //These are not needed in this function
  virtual void backup(FILE* file) {};
  virtual void restoreBackup(FILE* file) {};

  virtual const char* getName() const { return "FrequencyExcitation"; };

 private:
  double freq_GHz, omega_dimless;
  double UNz_ampl;
};

#endif
