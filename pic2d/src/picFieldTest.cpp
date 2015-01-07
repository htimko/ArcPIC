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

  picFieldTest.cpp:
  Test routine calculating the error from estimating the field at r=z=0
  using PIC instead of directly

***********************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

#include  <slu_ddefs.h>

#define   XTRN
#include  "pic.h"
#include  "var.h"
#include  "dim.h"
#include  "arrays1.h"
#include  "outp.h"
#include  "mydef.h"

//#include  "random.h"

#include  "init.h"
#include  "my_mem.h"
//#include  "my_time.h"
#include  "phi.h"
#include  "efield.h"
//#include  "push.h"
#include  "density.h"
//#include  "e_ion.h"
//#include  "moms.h"
//#include  "aver.h"
#include  "filenames.h"
#include  "outputz.h"
//#include  "engy.h"
//#include  "colls.h"
//#include  "order.h"
#include  "input.h"
//#include  "vdf.h"
//#include  "backup.h"
#include  "print_par.h"
//#include  "checkbounds.h"


int main(int argc, char* argv[]) {
  
  //input arguments
  if (not (argc == 5 || argc == 7)) {
    printf("Usage: picFieldTest mirror_order nzMax nrMax pointsPerGrid (zWitness=0.0 rWitness=0.0) \n");
    exit(0);
  }
  int mirror_order  = atoi(argv[1]);
  int nzMax         = atoi(argv[2]);
  int nrMax         = atoi(argv[3]);
  int pointsPerGrid = atoi(argv[4]);
  
  double zWitness(0.0), rWitness(0.0);
  if (argc == 7) {
    zWitness = atof(argv[5]);
    rWitness = atof(argv[6]);
  }
  
  printf("Got arguments:\n");
  printf(" - mirror_order  = %i\n", mirror_order);
  printf(" - nzMax         = %i\n", nrMax);
  printf(" - nrMax         = %i\n", nzMax);
  printf(" - pointsPerGrid = %i\n", pointsPerGrid);
  printf(" - zWitness      = %g\n", zWitness);
  printf(" - rWitness      = %g\n", rWitness);
  printf("\n");
  
  if (mirror_order < 0) {
    printf("Mirror order should be >= 0!\n");
    exit(1);
  }

  // SuperLU parameters
  SuperMatrix L_slu, U_slu;
  int *perm_c_slu;
  int *perm_r_slu;
  double *rhs_slu;

  int NG = NGR*NGZ;  

  //Read "input" (even if it loads some extra/unneeded stuff)
  input();

  double U0 = 0;
  double UNz = 0;

  //Sanity check of dim.h
  if (NGR != nr+1 || NGZ != nz+1) {
    printf("ERROR!! (NGR,NGZ) = (%i,%i), while (nr,nz)=(%i,%i)! (should be one smaller)\n", NGR, NGZ,nr,nz);
    printf("Aborting!\n");
    exit(1);
  }

  dr = dz;
  
  //Initialize
  allocate_arrays( nr, nz, &perm_c_slu, &perm_r_slu, &rhs_slu );
  calc_parameters_2D();
  M_ions[0] =  mi_over_me;                 //  H+
  M_ions[1] =  63.546/1.0073*mi_over_me-1; //  Cu+
  M_ions[2] =  63.546/1.0073*mi_over_me;   //  Cu
  
  for (int sort=0; sort < NSpecies; sort++) {
    cs_ions[sort] =  sqrt(M_ions[0]/M_ions[sort])*vi_0;
    vt_ions[sort] =  sqrt(M_ions[0]/M_ions[sort])*v_ti;
  }
  
  if (not (BC == 0 || BC == 1 || BC == 4)) {
    printf("Wrong BC.\n");
    exit(1);
  }
  
  potential_factorise_2D( nr, nz, NR, NZ, dr, dz, &L_slu, &U_slu, &perm_c_slu, &perm_r_slu );

  for (int i=0; i<NSpecies*NGR*NGZ; i++) {
    n_i[i]=0.;
  } 
  for (int i=0; i<Lastion*NGR*NGZ; i++) {
    E_ion_r[i]=0.;
    E_ion_z[i]=0.;
  } 
  for (int i=0; i<NGR*NGZ; i++) {
    n_e[i]=0.;
    E_grid_r[i]=0.;
    E_grid_z[i]=0.;
    phi[i]=0.;
  }
  
  // Initialise particles
  nr_e = 0, nr_i[0] = 0, nr_i[1] = 0, nr_i[2] = 0;

  init_reactions();
  print_parameters_2D();
  
  //Dummy, meant to keep analyses happy
  FILE* timeIndex;
  timeIndex = fopen("out/timeIndex.dat", "w");
  fprintf(timeIndex, "# StepNum SimTime[ns]\n");
  fflush(timeIndex);  
  
  //Get potential
  file_names_2D( 1 );
  //out_phi_2D(phi, 1, nr, nz, NZ, Omega_pe, dr, dz, fphi);
  //out_efield_2D( E_grid_z, E_grid_r, 1, nr, nz, NZ, Omega_pe, dr, dz, fEz, fEr );
  fprintf(timeIndex, "%08i %010e\n", 0, 0.0 );
  fflush(timeIndex);


  potential_backsolve_2D( nr, nz, NR, NZ, dz, U0, UNz, 
			  phi, L_slu, U_slu, perm_c_slu, perm_r_slu,
			  n_e, n_i + NG, &rhs_slu );
  electric_field_2D( phi, E_grid_r, E_grid_z, E_ion_r, E_ion_z, nr, nz, NR, NZ );
  density_2D( n_e, nr, nz, NR, NZ, Vcell, elec, qe, nr_e, 0, 1 );  

  file_names_2D( 2 );
  //out_phi_2D(phi, 1, nr, nz, NZ, Omega_pe, dr, dz, fphi);
  //out_efield_2D( E_grid_z, E_grid_r, 1, nr, nz, NZ, Omega_pe, dr, dz, fEz, fEr );
  fprintf(timeIndex, "%08i %010e\n", 1, 1.0 );
  fflush(timeIndex);

  double fieldConst = N_sp*1.602176565e-19/(4*PI*8.854187e-12); // N_sp*e/(4*pi*eps0) [V*m]
  double unitConst0 = dz*Ldb*1e-2;                         //dimless length (dz) -> m
  double unitConst1 = T_ref * SQU(dz/Omega_pe);            //dimless potential -> V
  double unitConst2 = unitConst1 * 2.0 / (Ldb*1e-2) / dz;  //dimless field -> V/m

  if (nrMax > nr) {
    printf("nrMax > nr, setting equal to nr = %i\n", nr);
    nrMax = nr;
  }
  if (nzMax > nz) {
    printf("nzMax > nz, setting equal to nz = %i\n", nz);
    nzMax = nz;
  }

  if (rWitness > nr) {
    printf("rWitness=%g > nr=%i, exiting\n", rWitness, nr);
    exit(1);
  }
  if (zWitness > nz) {
    printf("zWitness=%g > nz=%i, exiting\n", zWitness, nz);
    exit(1);
  }

  if (rWitness != 0.0) {
    printf("No analytical solution known for rWitness != 0.0, it will be incorrect -- now assuming rWitness=0 for analytics\n");
  }

  double witness_hr  = rWitness;
  int    witness_j   = (int)witness_hr;
  witness_hr        -= witness_j;       
  
  double witness_hz  = zWitness;
  int    witness_k   = (int)witness_hz;
  witness_hz        -= witness_k;    

  //Loop over positions to put test charges
  FILE* testOfile = fopen("tests/picFieldTest.dat", "w");
  fprintf(testOfile, "# idx iz ir z r fieldPIC[V/m] fieldDIRECT_z[V/m] fieldDIRECT_r[V/m]\n");
  fprintf(testOfile, "# nzMax=%i nrMax=%i pointsPerGrid=%i mirror_order=%i zWitness=%g rWitness=%g \n",
	  nzMax, nrMax, pointsPerGrid, mirror_order, zWitness, rWitness);

  int idx = 3;
  printf("Solving %i times on zXr = %i x %i...\n",
	 nzMax*nrMax*pointsPerGrid*pointsPerGrid,
	 nzMax*pointsPerGrid, nrMax*pointsPerGrid );
  for (int iz = 0; iz < nzMax*pointsPerGrid; iz++) {
    printf("iz = %i / %i \n", iz, nzMax*pointsPerGrid);
    for (int ir = 0; ir < nrMax*pointsPerGrid; ir++) {
      double z = iz/(1.0*pointsPerGrid);
      double r = ir/(1.0*pointsPerGrid);
      
      nr_e = 1;
      elec[0].p.z = z;
      elec[0].p.r = r;
      elec[0].p.vz = 0.0;
      elec[0].p.vr = 0.0;
      elec[0].p.vt = 0.0;
      elec[0].p.m = 1;
      
      if (z > nz || r > nr) {
	printf("Error: particle out of bounds! z = %f, r=%f\n", z, r);
	exit(1);
      }

      //Numerical solution
      density_2D( n_e, nr, nz, NR, NZ, Vcell, elec, qe, nr_e, 0, 1 );
      potential_backsolve_2D( nr, nz, NR, NZ, dz, U0, UNz, 
			      phi, L_slu, U_slu, perm_c_slu, perm_r_slu,
			      n_e, n_i + NG, &rhs_slu );
      electric_field_2D( phi, E_grid_r, E_grid_z, E_ion_r, E_ion_z, nr, nz, NR, NZ );

      double fieldPIC_z = (1-witness_hr)*(1-witness_hz)*E_grid_z[witness_j*NZ+witness_k] +
	                  (1-witness_hr)*witness_hz*E_grid_z[witness_j*NZ+(witness_k+1)] +
	                  witness_hr*(1-witness_hz)*E_grid_z[(witness_j+1)*NZ+witness_k] + 
	                  witness_hr*witness_hz*E_grid_z[(witness_j+1)*NZ+(witness_k+1)]; 
      double fieldPIC_r = (1-witness_hr)*(1-witness_hz)*E_grid_r[witness_j*NZ+witness_k] + 
	                  (1-witness_hr)*witness_hz*E_grid_r[witness_j*NZ+(witness_k+1)] + 
	                  witness_hr*(1-witness_hz)*E_grid_r[(witness_j+1)*NZ+witness_k] + 
	                  witness_hr*witness_hz*E_grid_r[(witness_j+1)*NZ+(witness_k+1)];

      //Field output at this position (uses LOTS of disk space -- easily 100's of GB -> slow)
      file_names_2D( idx );
      //out_phi_2D(phi, 1, nr, nz, NZ, Omega_pe, dr, dz, fphi);
      //out_efield_2D( E_grid_z, E_grid_r, 1, nr, nz, NZ, Omega_pe, dr, dz, fEz, fEr );
      //out_dens_2D( n_e, 1,    -1., nr, nz, NZ, Omega_pe, dr, dz, fn_e );      
      //Index file needed for standard analysis
      fprintf(timeIndex, "%08i %010e\n", idx-1, idx-1.0 );
      fflush(timeIndex);
      
      //Analytical solution
      //Hold each term separately, then sum backwards
      double* fieldDirect_orders = new double[mirror_order+1]; 
      
      //0-order
      //The charge itself (negative => positive field)
      double R3   = pow( SQU(z-zWitness) + SQU(r), 3.0/2.0);
      fieldDirect_orders[0] = -1.0*fieldConst * ( -(z-zWitness) / (R3 * SQU(unitConst0)) );
      //Neuman boundary (r=nr) (INCORRECT!)
      //double r1 = 2*nr-r;
      //R3   = pow( SQU(z) + SQU(r1), 3.0/2.0);
      //fieldDirect += -1.0*fieldConst * ( -z / (R3 * SQU(unitConst0)) );
      
      //Higher orders dirichlet boundary
      for (int o = 1; o <= mirror_order; o++) {
	int chargeSign = ((o % 2)*2 - 1); //+1 on odd indexes, -1 on even indexes
	
	//cathode side mirror charges
	double z1 = -o*nz + ( o%2==1 ? nz-z : z);
	//printf("o=%i chargeSign=%i z1=%f ", o, chargeSign, z1);
	R3   = pow( SQU(z1-zWitness) + SQU(r),  3.0/2.0);
	fieldDirect_orders[o] = chargeSign*fieldConst * ( -(z1-zWitness) / (R3 * SQU(unitConst0)) );
	//Neumann boundary (r=nr) (INCORRECT!)
	//R3   = pow( SQU(z1) + SQU(r1), 3.0/2.0);
	//fieldDirect += chargeSign*fieldConst * ( -z1 / (R3 * SQU(unitConst0)) );
	
	//cathode side
	double z2 =  o*nz + ( o%2==1 ? nz-z : z);
	//printf("z2=%f\n",z2);
	R3   = pow( SQU(z2-zWitness) + SQU(r),  3.0/2.0);
	fieldDirect_orders[o] += chargeSign*fieldConst * ( -(z2-zWitness) / (R3 * SQU(unitConst0)) );
	//Neumann boundary (r=nr) (INCORRECT!)
	//R3   = pow( SQU(z2) + SQU(r1), 3.0/2.0);
	//fieldDirect += chargeSign*fieldConst * ( -z2 / (R3 * SQU(unitConst0)) );

      }

      //Sum up the total field - backwards!
      for (int o = mirror_order; o >= 0; o--) {
	for (int o2 = o-1; o2 >= 0; o2--) {
	  fieldDirect_orders[o] += fieldDirect_orders[o2];
	}
      }
      double fieldDirect = fieldDirect_orders[mirror_order];

      fprintf(testOfile, "%i %i %i %f %f %f %f %f",
	      idx, iz, ir, z, r, fieldPIC_z*unitConst2, fieldPIC_r*unitConst2, fieldDirect);
      for (int o = 0; o <= mirror_order; o++) {
	fprintf(testOfile, " %f", fieldDirect_orders[o]);
      }
      fprintf(testOfile, "\n");
      fflush(testOfile);

      delete[] fieldDirect_orders;

      idx++;
    }
  }

  fclose(testOfile);
  fclose(timeIndex);
  
}
