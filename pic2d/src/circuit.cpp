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

  circuit.cpp:
  Models for the external circuit.

***********************************************************************/

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "circuit.h"

#define XTRN extern
#include "mydef.h"
#include "pic.h" //Needed for var.h
#include "var.h"
#undef XTRN

#include "input.h"

using namespace std;

// ******** Implementation of Circuit ******************
Circuit* Circuit::LoadCircuit(FILE* in_file) {
  vector<char*> options;
  char* circuitName = readInputSection(in_file, options);
  
  Circuit* ret = NULL;

  if ( ! strncmp(circuitName, "FixedVoltage", NAME_MAXLEN) ) {
    ret = new FixedVoltage(options);
  }
  else if ( ! strncmp(circuitName, "FixedVoltage_resistor", NAME_MAXLEN) ) {
    ret = new FixedVoltage_resistor(options);
  }
  else if ( ! strncmp(circuitName, "FixedVoltage_resistorCapacitor", NAME_MAXLEN) ) {
    ret = new FixedVoltage_resistorCapacitor(options);
  }
  else if ( ! strncmp(circuitName, "FixedVoltage_ramp", NAME_MAXLEN) ) {
    ret = new FixedVoltage_ramp(options);
  }
  else if ( ! strncmp(circuitName, "TwoCapacitors", NAME_MAXLEN) ) {
    ret = new TwoCapacitors(options);
  }
  else if ( ! strncmp(circuitName, "FrequencyExcitation", NAME_MAXLEN) ) {
    ret = new FrequencyExcitation(options);
  }
  else {
    cout << "Error in LoadCircuit: Can't handle given circuitName '"
	 << circuitName << "'" << endl;
    exit(1);
  }
  
  //Cleanup
  delete[] circuitName;
  for (size_t i = 0; i < options.size(); i++) {
    delete[] options[i];
  }
  options.clear();
  
  return ret;
}

void Circuit::initFile(){
  ofile = fopen("circuit.dat", "w");
  fprintf(ofile, "## Circuit : %s ##\n", this->getName());
  fprintf(ofile, "## timestep deltaQ_accumulated[#superparticles] deltaU[dimless]\n");
  fflush(ofile);
}
void Circuit::writeFile(unsigned int nstep){
  if (ofile == NULL) initFile();
  if ( nstep % file_timestep == 0) {
    fprintf(ofile, "%u %10f %10f\n", nstep, deltaQ_accumulated, getUNz()-getU0());
    deltaQ_accumulated = 0;
    fflush(ofile);
  }
}

// ******** Implementation of FixedVoltage ******************

FixedVoltage::FixedVoltage( std::vector<char*> options ) { 
  if (options.size() != 3) {
    cout << "Error in FixedVoltage(): Expected 3 options, got "
	 << options.size() << endl;
    exit(1);
  }

  //Interpret lines
  sscanf(options[0], "%*[^:]%*[:] %lg", &(this->U0));
  sscanf(options[1], "%*[^:]%*[:] %lg", &(this->UNz));
  sscanf(options[2], "%*[^:]%*[:] %u",  &(this->file_timestep) );
}

void FixedVoltage::init() {
  this->U0  *= SQU(Omega_pe/dz);
  this->UNz *= SQU(Omega_pe/dz);
}

void FixedVoltage::print_par() {
  // Need to convert back from dimless
  printf( " - U0:                  %g [T_ref]\n",
	  (this->U0)*dz*dz/Omega_pe/Omega_pe );
  printf( " - UNz:                 %g [T_ref]\n",
	  (this->UNz)*dz*dz/Omega_pe/Omega_pe );
  printf( " - file_timestep:       %u\n",
	  this->file_timestep );
}

void FixedVoltage::timestep(double deltaQ, unsigned int nstep) {
  //Write to file
  deltaQ_accumulated += deltaQ;
  writeFile(nstep);
}

// ******** Implementation of FixedVoltage_resistor ******************

FixedVoltage_resistor::FixedVoltage_resistor( std::vector<char*> options ) { 
  if (options.size() != 3) {
    cout << "Error in FixedVoltage_resistor(): Expected 3 options, got "
	 << options.size() << endl;
    exit(1);
  }

  this->U0 = 0.0;
  //Interpret lines
  sscanf(options[0], "%*[^:]%*[:] %lg", &(this->UNz));
  this->Vin = this->UNz - this->U0;
  sscanf(options[1], "%*[^:]%*[:] %lg", &(this->Rseries_dim));
  sscanf(options[2], "%*[^:]%*[:] %u",  &(this->file_timestep) );
}

void FixedVoltage_resistor::init() {
  this->U0      *= SQU(Omega_pe/dz);
  this->UNz     *= SQU(Omega_pe/dz);
  this->Vin     *= SQU(Omega_pe/dz);
  this->Rseries  = this->Rseries_dim * 3.713294e-6*Omega_pe*sqrt(T_ref)/SQU(dz)/Ndb;
}

void FixedVoltage_resistor::print_par() {
  // Need to convert back from dimless
  printf( " - U0:                  %g [T_ref]\n",
	  (this->U0)*dz*dz/Omega_pe/Omega_pe );
  printf( " - UNz:                 %g [T_ref]\n",
	  (this->UNz)*dz*dz/Omega_pe/Omega_pe );
  printf( " - Vin:                 %g [T_ref]\n",
	  (this->Vin)*dz*dz/Omega_pe/Omega_pe );
  printf( " - Rseries              %g [Ohm] = %g [1] \n",
	  this->Rseries_dim, this->Rseries);
  printf( " => Assuming burning voltage = 50 V, max current = %g A\n",
	  (Vin*SQU(dz/Omega_pe)*T_ref - 50.0)/Rseries_dim );
  printf( " - file_timestep:       %u\n",
	  this->file_timestep );
}

void FixedVoltage_resistor::timestep(double deltaQ, unsigned int nstep) {
  UNz = Vin - deltaQ*Rseries;
  
  //Write to file
  deltaQ_accumulated += deltaQ;
  writeFile(nstep);
}

// ******** Implementation of FixedVoltage_resistorCapacitor ******************

FixedVoltage_resistorCapacitor::FixedVoltage_resistorCapacitor( std::vector<char*> options ) { 
  if (options.size() != 4) {
    cout << "Error in FixedVoltage_resistorCapacitor(): Expected 4 options, got "
	 << options.size() << endl;
    exit(1);
  }

  this->U0 = 0.0;
  
  //Interpret lines
  sscanf(options[0], "%*[^:]%*[:] %lg", &(this->UNz));
  this->Vin = this->UNz - this->U0;
  sscanf(options[1], "%*[^:]%*[:] %lg", &(this->Rseries_dim));
  sscanf(options[2], "%*[^:]%*[:] %lg", &(this->Cgap_dim));
  sscanf(options[3], "%*[^:]%*[:] %u",  &(this->file_timestep) );
}

void FixedVoltage_resistorCapacitor::init() {
  this->U0      *= SQU(Omega_pe/dz);
  this->UNz     *= SQU(Omega_pe/dz);
  this->Vin     *= SQU(Omega_pe/dz);
  this->Rseries  = this->Rseries_dim * 3.713294e-6*Omega_pe*sqrt(T_ref)/SQU(dz)/Ndb;
  
  this->Cgap     = this->Cgap_dim    * 1.519260e10*SQU(dz/Omega_pe)*Ndb*sqrt(n_ref/T_ref);
  
  this->I_circ   = 0.0;
}

void FixedVoltage_resistorCapacitor::print_par() {
  // Need to convert back from dimless
  printf( " - U0:                  %g [T_ref]\n",
	  (this->U0)*dz*dz/Omega_pe/Omega_pe );
  printf( " - UNz:                 %g [T_ref]\n",
	  (this->UNz)*dz*dz/Omega_pe/Omega_pe );
  printf( " - Vin:                 %g [T_ref]\n",
	  (this->Vin)*dz*dz/Omega_pe/Omega_pe );
  printf( " - Rseries:             %g [Ohm] = %g [1] \n",
	  this->Rseries_dim, this->Rseries);
  printf( " => Assuming burning voltage = 50 V, max current = %g A\n",
	  (Vin*SQU(dz/Omega_pe)*T_ref - 50.0)/Rseries_dim );
  printf( " - Cgap:                %g [Farad]\n",
	  this->Cgap );
  printf( " - file_timestep:       %u\n",
	  this->file_timestep );
}

void FixedVoltage_resistorCapacitor::timestep(double deltaQ, unsigned int nstep) {
  UNz = UNz + (I_circ-deltaQ)/Cgap;
  I_circ = (Vin-UNz)/Rseries;
  
  //Write to file
  deltaQ_accumulated += deltaQ;
  writeFile(nstep);
}

void FixedVoltage_resistorCapacitor::initFile(){
  ofile = fopen("circuit.dat", "w");
  fprintf(ofile, "## Circuit : %s ##\n", this->getName());
  fprintf(ofile, "## timestep deltaQ_accumulated[#superparticles] deltaU[dimless] | I_circ[superparticles/timestep] \n");
  fflush(ofile);
}
void FixedVoltage_resistorCapacitor::writeFile(unsigned int nstep){
  if (ofile == NULL) initFile();
  if ( nstep % file_timestep == 0) {
    fprintf(ofile, "%u %10f %10f %10f\n", nstep, deltaQ_accumulated, getUNz()-getU0(), I_circ);
    deltaQ_accumulated = 0;
    fflush(ofile);
  }
}

// ******** Implementation of FixedVoltage_ramp ******************

FixedVoltage_ramp::FixedVoltage_ramp( std::vector<char*> options ) {
  if (options.size() != 4) {
    cout << "Error in FixedVoltage_ramp(): Expected 4 options, got "
	 << options.size() << endl;
    exit(1);
  }

  //Interpret lines
  sscanf(options[0], "%*[^:]%*[:] %lg", &(this->U0));
  sscanf(options[1], "%*[^:]%*[:] %lg", &(this->UNz));
  sscanf(options[2], "%*[^:]%*[:] %lg", &(this->rampTime));
  sscanf(options[3], "%*[^:]%*[:] %u",  &(this->file_timestep) );
}
void FixedVoltage_ramp::init() {
  this->U0  *= SQU(Omega_pe/dz);
  this->UNz *= SQU(Omega_pe/dz);
  this->Vfinal = this->UNz - this->U0;
  this->UNz = U0;
  this->rampTime_steps = 
    ceil( 56414.6*sqrt(n_ref)/Omega_pe * 1e-9*(this->rampTime) );
}
void FixedVoltage_ramp::timestep(double deltaQ, unsigned int nstep) {
  if (nstep < this->rampTime_steps) {
    this->UNz = ( double(nstep)/double(this->rampTime_steps) ) * this->Vfinal + this->U0;
  }
  else {
    this->UNz = U0 + Vfinal;
  }

  //Write to file
  deltaQ_accumulated += deltaQ;
  writeFile(nstep);
}

void FixedVoltage_ramp::print_par() {
  // Need to convert back from dimless
  printf( " - U0:                  %g [T_ref]\n",
	  (this->U0)*dz*dz/Omega_pe/Omega_pe );
  printf( " - UNz:                 %g [in T_ref]\n",
	  (this->UNz)*dz*dz/Omega_pe/Omega_pe );
  printf( " - rampTime             %g [ns] = %u [steps]\n",
	  this->rampTime, this->rampTime_steps );
  printf( " - file_timestep:       %u\n",
	  this->file_timestep );
}


// ************* Implementaiton of TwoCapacitors ******************
TwoCapacitors::TwoCapacitors( std::vector<char*> options ) {
  if (options.size() != 8) {
    cout << "Error in TwoCapacitors(): Expected 8 options, got "
	 << options.size() << endl;
    exit(1);
  }
  
  //Interpret lines
  // (save in [dim] variables, convert to dimless later.
  // These variables are NOT UPDATED! )
  sscanf(options[0], "%*[^:]%*[:] %lg", &(this->U0_dim));
  sscanf(options[1], "%*[^:]%*[:] %lg", &(this->UNz_dim));
  sscanf(options[2], "%*[^:]%*[:] %lg", &(this->C_ext_dim));
  sscanf(options[3], "%*[^:]%*[:] %lg", &(this->C_ext2_dim));
  sscanf(options[4], "%*[^:]%*[:] %lg", &(this->R_ext_dim));
  
  sscanf(options[5], "%*[^:]%*[:] %lg", &(this->t_change_dim));
  if (this->t_change_dim < 0.0) {
    cout << "Error in TwoCapacitors()::init: t_change must be > 0.0" << endl;
  }
  
  char foo;
  sscanf(options[6],  "%*[^:]%*[:] %c",  &(foo));
  if      (foo == 'y') this->avoid_negative = true;
  else if (foo == 'n') this->avoid_negative = false;
  else {
    printf("Error in TwoCapacitors::TwoCapacitors(): avoid_negative has to be either 'y' or 'n' \n");
    exit(1);
  }  

  
  sscanf(options[7], "%*[^:]%*[:] %u",  &(this->file_timestep) );
}

void TwoCapacitors::init() {
  this->U0     = U0_dim  * SQU(Omega_pe/dz);
  this->UNz    = UNz_dim * SQU(Omega_pe/dz);
  
  this->C_ext  = C_ext_dim  * 1.519260e10*SQU(dz/Omega_pe)*Ndb*sqrt(n_ref/T_ref);
  this->C_ext2 = C_ext2_dim * 1.519260e10*SQU(dz/Omega_pe)*Ndb*sqrt(n_ref/T_ref);
  this->R_ext  = R_ext_dim  * 3.713294e-6*Omega_pe*sqrt(T_ref)/SQU(dz)/Ndb;
  this->Q_ext  = (this->C_ext)*(UNz-U0);

  this->t_change = (unsigned int) ( t_change_dim*1e-9*56414.6*sqrt(n_ref) / Omega_pe );

  this->check_C = false;
}
void TwoCapacitors::re_init() {
  ofile = fopen("circuit.dat", "a");
  
  C_ext  = C_ext_dim  * 1.519260e10*SQU(dz/Omega_pe)*Ndb*sqrt(n_ref/T_ref);
  C_ext2 = C_ext2_dim * 1.519260e10*SQU(dz/Omega_pe)*Ndb*sqrt(n_ref/T_ref);
  R_ext  = R_ext_dim  * 3.713294e-6*Omega_pe*sqrt(T_ref)/SQU(dz)/Ndb;
  Q_ext *=  1.519260e10*Ndb*sqrt(n_ref/T_ref)/T_ref;

  this->t_change = (unsigned int) ( t_change_dim*1e-9*56414.6*sqrt(n_ref) / Omega_pe );

  //TODO: Uses global variable (disable unused code for now)
  //if ( ( nsteps > t_change ) ) {
  //  check_C = false;
  //}
}

void TwoCapacitors::print_par() {
  printf( " - U0:                  %g [T_ref]\n",  this->U0_dim );
  printf( " - UNz:                 %g [T_ref]\n",  this->UNz_dim);
  printf( " - C_ext:               %g [Farad]\n",  this->C_ext_dim );
  printf( " - C_ext2:              %g [Farad]\n",  this->C_ext2_dim );
  printf( " - R_ext:               %g [Ohm]\n",    this->R_ext_dim );
  printf( " - t_change             %g [ns] = %u [dt]\n", this->t_change_dim, this->t_change);
  printf( " - avoid_negative       %c \n",         this->avoid_negative ? 'y' : 'n' );
  printf( " - file_timestep:       %u\n",          this->file_timestep );
}

void TwoCapacitors::timestep(double deltaQ, unsigned int nstep) {
  // after 5 ns, use external capacitor capacitance
  if ( (nstep > t_change) && (this->check_C == false) ) { 
    this->UNz   = (this->Q_ext)/(this->C_ext) + (this->U0);
    this->C_ext = this->C_ext2;
    this->Q_ext = this->C_ext*(this->UNz - this->U0);
    check_C = true;
  }

  //Update system state
  this->Q_ext -= deltaQ;
  this->UNz = (this->Q_ext)/(this->C_ext) - (this->R_ext)*deltaQ + this->U0;
  
  //Don't allow negative voltage accross the gap
  if ( avoid_negative and this->UNz < this->U0 ) this->UNz = this->U0;

  //Write to file
  deltaQ_accumulated += deltaQ;
  writeFile(nstep);
}

void TwoCapacitors::backup(FILE* file) {
  // In SI-units 
  fprintf(file,"%19.11e %19.11e %19.11e %19.11e\n",
	  (this->R_ext)  / 3.713294e-6/Omega_pe/sqrt(T_ref)*SQU(dz)*Ndb, 
	  (this->C_ext)  / 1.519260e10/SQU(dz/Omega_pe)/Ndb/sqrt(n_ref/T_ref),
	  (this->C_ext2) / 1.519260e10/SQU(dz/Omega_pe)/Ndb/sqrt(n_ref/T_ref),
	  (this->Q_ext)  / 1.519260e10/Ndb/sqrt(n_ref/T_ref)*T_ref );
}
void TwoCapacitors::restoreBackup(FILE* file) {
  fscanf(file,"%19lg %19lg %19lg %19lg",
	 &(this->R_ext_dim), &(this->C_ext_dim), &(this->C_ext2_dim), &(this->Q_ext) );
}

void TwoCapacitors::initFile(){
  ofile = fopen("circuit.dat", "w");
  fprintf(ofile, "## Circuit : %s ##\n", this->getName());
  fprintf(ofile, "## timestep deltaQ_accumulated[#superparticles] deltaU[dimless] | Q_ext check_C\n");
  fflush(ofile);
}
void TwoCapacitors::writeFile(unsigned int nstep){
  if (ofile == NULL) initFile();
  if ( nstep % file_timestep == 0) {
    fprintf(ofile, "%u %10f %10f %10f %c \n", nstep, deltaQ_accumulated, getUNz()-getU0(), Q_ext, check_C ? 'y' : 'n');
    deltaQ_accumulated = 0;
    fflush(ofile);
  }
}

// ******** Implementation of FrequencyExcitation ******************

FrequencyExcitation::FrequencyExcitation( std::vector<char*> options ) { 
  if (options.size() != 3) {
    cout << "Error in FrequencyExcitation(): Expected 3 options, got "
	 << options.size() << endl;
    exit(1);
  }

  //Interpret lines
  sscanf(options[0], "%*[^:]%*[:] %lg", &(this->freq_GHz));
  sscanf(options[1], "%*[^:]%*[:] %lg", &(this->UNz_ampl));
  sscanf(options[2], "%*[^:]%*[:] %u",  &(this->file_timestep) );
}

void FrequencyExcitation::init() {
  this->U0             = 0.0;
  this->UNz            = 0.0;
  this->UNz_ampl      *= SQU(Omega_pe/dz);
  this->omega_dimless  = TWOPI*freq_GHz*1e9 * Omega_pe/(56313.6*sqrt(n_ref));
}

void FrequencyExcitation::print_par() {
  // Need to convert back from dimless
  printf( " - UNz_ampl:            %g [T_ref]\n",
	  (this->UNz_ampl)*dz*dz/Omega_pe/Omega_pe );
  printf( " - freq_GHz =           %g [GHz] = %g [radians/timestep]\n", 
	  freq_GHz, omega_dimless);
  printf( " - file_timestep:       %u\n",
	  this->file_timestep );
}

void FrequencyExcitation::timestep(double deltaQ, unsigned int nstep) {
  
  this->UNz = UNz_ampl*sin(omega_dimless*nstep);
  
  //Write to file
  deltaQ_accumulated += deltaQ;
  writeFile(nstep);
}
