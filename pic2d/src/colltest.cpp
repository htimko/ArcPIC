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

  colltest.cpp:
  Test routine for checking that collisions work properly.

***********************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#include "random.h";

#define   XTRN
#include  "pic.h"
#include  "var.h"
#include "colls.h"

int main(int argc, char* argv[]) {
  //Initialize scaling etc.
  T_ref = 4.0; //eV
  n_ref = 1e10; //cm^-3
  Ndb = n_ref*pow(7.43e2*sqrt(T_ref/n_ref),3)/4.0; //Ldb^-3 //Want N_sp=4

  seed = 1234;
  dt_ion = 1;
  Omega_pe = 0.2;
  dz = 0.5;
  Ldb = 7.43e2*sqrt(T_ref/n_ref);
  N_sp = n_ref*Ldb*Ldb*Ldb/Ndb;
  ncoll_el = 1; //Needed for OLD version
  dt_ion = 5; //This will show up if its there

  double vte = sqrt(T_ref/511e3)*3e8; //Thermal velocity
  double vk = sqrt(2)*vte; //Mode of P(|v|)
  double Emean = 3.0*T_ref/2.0; //Mean energy
  cout << "vte = " << vte << ", vk = " << vk << " [m/s]" << endl;
  
  double vconv = Omega_pe/(dz*vte); //v -> \tilde v
  cout << "vconv = " << vconv << endl;

  //cout << "vk -> eV: " << 0.5*vk*vk*511e3/3e8/3e8 << endl; //OK
  
  int numElec = n_ref*M_PI*pow(0.5*sqrt(552635*T_ref/n_ref),3)/N_sp;
  unsigned int numSteps = 4*141980; // 4*nu0^-1

  cout << "N_sp = " << N_sp << ", numElec = " << numElec << endl;

  //System initialization
  rng_initialize(seed, 1);
  Particle* elec = new Particle[numElec];
  int e_order[1] = {numElec};
  Vec3d mcheck = {0.0,0.0,0.0};        
  double echeck = 0.0;
  
  //Initialize velocity distribution to uniform
  // energy distribution with mean E such that
  // it corresponds to the correct temperature
  for (int i = 0; i < numElec; i++) {
    double E = Emean*(1+0.2*2*(RAND-0.5));
    double v = sqrt(2*E/511e3)*3e8*vconv;
    
    //Direction: uniform distribution over sphere
    double theta = acos(2*RAND-1);
    double phi   = RAND*2*M_PI;
    elec[i].p.vz = v*cos(theta);
    elec[i].p.vr = v*sin(theta)*cos(phi);
    elec[i].p.vt = v*sin(theta)*sin(phi);

    //cout << sqrt(elec[i].p.vz*elec[i].p.vz + elec[i].p.vr*elec[i].p.vr + elec[i].p.vt*elec[i].p.vt)/(vk*vconv) << endl; //ok
  }

  //Start collisions
  cout << "Colliding " << numElec << " electrons..." << endl;
  ofstream eOfile;
  eOfile.open("energy.out");
  eOfile.precision(16);
  for (unsigned int i = 0; i < numSteps; i++) {
    //Particle velocities for calculating the velocity distribution
    if (i % 5000 == 0) {
      
      ofstream ofile;
      char ofname[50];
      sprintf(ofname, "vdist_%09u.out", i);
      cout << "\t Output: " << ofname << endl;
      ofile.open(ofname);
      for (int j = 0; j < numElec; j++) {
	ofile << elec[j].p.vz/vconv << ", "
	      << elec[j].p.vr/vconv << ", "
	      << elec[j].p.vt/vconv << endl;
      }
      ofile.close();
    }

    //Total energy [eV]
    if (i % 100 == 0) {
      double E = 0.0;
      for (int j = 0; j < numElec; j++) {
	E += elec[j].p.vz*elec[j].p.vz/vconv/vconv +
	  elec[j].p.vr*elec[j].p.vr/vconv/vconv +
	  elec[j].p.vt*elec[j].p.vt/vconv/vconv;
	/*
	  cout << sqrt(elec[j].p.vz*elec[j].p.vz/vconv/vconv +
	  elec[j].p.vr*elec[j].p.vr/vconv/vconv +
	  elec[j].p.vt*elec[j].p.vt/vconv/vconv) / vk << endl;
	*/
      }
      E *= 0.5*511e3*E/3e8/3e8;
      eOfile << i << " " << E << " " << E / ((double)numElec) << endl;

      cout << "step = " << i << endl;
      cout << "\t echeck = sum_steps (ebefore-eafter) = " << echeck << endl;
      cout << "\t momcheck (z,r,t) = (" 
	   << mcheck.z << ", " << mcheck.r << ", " << mcheck.t << ")" << endl;
    }  

    coll_el_knm_2D(elec, numElec, e_order, 1, 1, 2, 1, 0, &mcheck, &echeck, 1);
    //coll_el_knm_2D_OLD(elec, numElec, e_order, 1, 1, 2, 1, 0, &mcheck, &echeck);
        
  }
  
  eOfile.close();
  
}
