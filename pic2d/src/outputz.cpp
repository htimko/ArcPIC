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

  outputz.cpp:
  Output and diagnostic routines

***********************************************************************/

#include <stdio.h>

#include "dim.h"
#include "pic.h"
#include "mydef.h"



/***************************
 Outputting routines 
***************************/

void out_fv_along_2D( double fvz[], double fvr[], double fvabs[], int nr, int nz, 
		      const char *dat1, const char *dat2, const char *dat3 ) {
  FILE *file1, *file2, *file3;
  int j, k, n;
  int halfnr = nr/2;
  int halfnz = nz/2;
  
  file1 = fopen(dat1, "w");
  file2 = fopen(dat2, "w");
  file3 = fopen(dat3, "w");


  for (j=0; j<halfnr; j++) {
    for (k=0; k<halfnz; k++) {
      
      for (n=0; n < Nvdst; n++) {
	fprintf(file1, "% .5e ", fvz[n + Nvdst*(j*(nz/2)+k)]);
	fprintf(file2, "% .5e ", fvr[n + Nvdst*(j*(nz/2)+k)]);
	fprintf(file3, "% .5e ", fvabs[n + Nvdst*(j*(nz/2)+k)]);
      }
      
      fprintf(file1, "\n");
      fprintf(file2, "\n");
      fprintf(file3, "\n");
    }
  }
  
  fclose(file1);
  fclose(file2);
  fclose(file3);
  
}




// temporary for density
void out_tempn_2D( double n_el[], double n_ion[], int n_aver_e, int n_aver_ion, 
	       int nr, int nz, int NZ, double Vcell[], double Omega_pe, double qp, 
	       double dr, double dz, const char *dat_e, const char *dat_i ) {
  
  FILE *file_e, *file_i;
  double fi, fe;
  int j, k;  
  
  if (n_aver_e < 1)  n_aver_e = 1;
  if (n_aver_ion < 1)  n_aver_ion = 1;
  
  //fi = qp*SQU(1./Omega_pe)*TWOPI/(dz*n_aver_ion);   
  //fe = qp*SQU(1./Omega_pe)*TWOPI/(dz*n_aver_e); 
  fi = qp*SQU(1./Omega_pe)/n_aver_ion;   //CORR! 31.03.2010
  fe = qp*SQU(1./Omega_pe)/n_aver_e;   //CORR! 31.03.2010  

  file_e = fopen(dat_e, "w");
  file_i = fopen(dat_i, "w");

  for (j=0; j<nr; j++) {
    for (k=0; k<nz; k++) {
      fprintf(file_e, "% .5e % .5e % .5e \n", dr*j, dz*k, fe*n_el[j*NZ+k]/Vcell[j]); 
      fprintf(file_i, "% .5e % .5e % .5e \n", dr*j, dz*k, fi*n_ion[j*NZ+k]/Vcell[j]); 
    }
  }
  
  fclose(file_e);
  fclose(file_i);
  
}

void out_n_2D( Moments mom_el[], Moments mom_ion[], int n_aver_e, int n_aver_ion, 
	       int nr, int nz, int NZ, double Vcell[], double Omega_pe, double qp, 
	       double dr, double dz, const char *dat_e, const char *dat_i ) {

  FILE *file_e, *file_i;
  double fi, fe;
  int j, k;  

  if (n_aver_e < 1)  n_aver_e = 1;
  if (n_aver_ion < 1)  n_aver_ion = 1;
  
  //fi = qp*SQU(1./Omega_pe)*TWOPI/(dz*n_aver_ion);   
  //fe = qp*SQU(1./Omega_pe)*TWOPI/(dz*n_aver_e);   
  fi = qp*SQU(1./Omega_pe)/n_aver_ion;   //CORR! 31.03.2010
  fe = qp*SQU(1./Omega_pe)/n_aver_e;   //CORR! 31.03.2010

  file_e = fopen(dat_e, "w");
  file_i = fopen(dat_i, "w");

  for (j=0; j<nr; j++) {
    for (k=0; k<nz; k++) {
      fprintf(file_e, "% .5e % .5e % .5e \n", dr*j, dz*k, fe*mom_el[j*NZ+k].n/Vcell[j]); 
      fprintf(file_i, "% .5e % .5e % .5e \n", dr*j, dz*k, fi*mom_ion[j*NZ+k].n/Vcell[j]); 
    }
  }
  
  fclose(file_e);
  fclose(file_i);
  
}




void out_dens_2D( double dens_av[], int n_aver, double sign, int nr, int nz, int NZ, double Omega_pe,  
		  double dr, double dz, const char *dat_p ) {
  FILE *file_p;
  double fp;
  int j, k;  
  
  if (n_aver < 1)  n_aver = 1;

  fp = sign*SQU(1./Omega_pe)/n_aver;  

  file_p = fopen(dat_p, "w");

  for (j=0; j<=nr; j++) {
    for (k=0; k<=nz; k++) {
      fprintf(file_p, "% .5e % .5e % .5e \n", dr*j, dz*k, fp*dens_av[j*NZ+k]); 
    }
  }
  
  fclose(file_p);
  
}

void out_phi_2D( double phi[], int n_aver, int nr, int nz, int NZ, 
		 double Omega_pe, double dr, double dz, const char *dat1 ) {

  FILE *file1;
  double f1;
  int j, k;
  
  if( n_aver < 1 )  n_aver = 1;
  
  f1 = SQU(dz/Omega_pe)/n_aver;   

  file1 = fopen(dat1, "w");
  for (j=0; j<=nr; j++) {
    for (k=0; k<=nz; k++) {
      fprintf(file1, "% .5e % .5e % .5e \n", dr*j, dz*k, f1*phi[j*NZ+k]); 
    }
  }
  fclose(file1); 
  
}




void out_efield_2D( double efield_z[], double efield_r[], int n_aver, int nr, int nz, int NZ, 
		    double Omega_pe, double dr, double dz, const char *dat1, const char *dat2 ) {

  FILE *file1, *file2;
  double f1;
  int j, k;
  
  if( n_aver < 1 )  n_aver = 1;
  
  f1 = 2*dz/SQU(Omega_pe)/n_aver;   

  file1 = fopen(dat1, "w");
  file2 = fopen(dat2, "w");

  for (j=0; j<=nr; j++) {
    for (k=0; k<=nz; k++) {
      fprintf(file1, "% .5e % .5e % .5e \n", dr*j, dz*k, f1*efield_z[j*NZ+k]); 
      fprintf(file2, "% .5e % .5e % .5e \n", dr*j, dz*k, f1*efield_r[j*NZ+k]); 
    }
  }
  fclose(file1); 
  fclose(file2);
  
}




void out_vels_2D( Moments  mom[], int nr, int nz, int NZ, double u0, 
		  double dr, double dz, const char *dat1, const char *dat2, const char *dat3 ) {

  FILE *file1, *file2, *file3;
  int j, k;
  
  u0 = MAX(u0,1.e-10);

  for (j=0; j<nr; j++) {
    for (k=0; k<nz; k++) {
      
      if( mom[j*NZ+k].n > 1 ) {
	mom[j*NZ+k].uz /= (u0*mom[j*NZ+k].n);
	mom[j*NZ+k].ur /= (u0*mom[j*NZ+k].n);
	mom[j*NZ+k].ut /= (u0*mom[j*NZ+k].n);
      }
      else {
	mom[j*NZ+k].uz = mom[j*NZ+k].ur = mom[j*NZ+k].ut = 0.;
      }      
    }
  }
  
  file1 = fopen(dat1, "w");
  file2 = fopen(dat2, "w");
  file3 = fopen(dat3, "w");
        
  for (j=0; j<nr; j++) {
    for (k=0; k<nz; k++) {
      fprintf(file1, "% .5e % .5e % .5e \n", dr*j, dz*k, mom[j*NZ+k].uz); 
      fprintf(file2, "% .5e % .5e % .5e \n", dr*j, dz*k, mom[j*NZ+k].ur);  
      fprintf(file3, "% .5e % .5e % .5e \n", dr*j, dz*k, mom[j*NZ+k].ut); 
    }
  }
  
  fclose(file1);
  fclose(file2);
  fclose(file3);
     	
}




void out_temps_2D( Moments  mom[], double u0, double fnorm, int nr, int nz, int NZ, 
		   double dr, double dz, const char *dat1, const char *dat2, const char *dat3 ) {

  FILE *file1, *file2, *file3;
  int j, k;
  
  for (j=0; j<nr; j++) {
    for (k=0; k<nz; k++) {
      
      if( mom[j*NZ+k].n > 1 ) {
	mom[j*NZ+k].tz = ( mom[j*NZ+k].tz/(u0*u0*mom[j*NZ+k].n) - SQU(mom[j*NZ+k].uz) )*fnorm;
	mom[j*NZ+k].tr = ( mom[j*NZ+k].tr/(u0*u0*mom[j*NZ+k].n) - SQU(mom[j*NZ+k].ur) )*fnorm;
	mom[j*NZ+k].tt = ( mom[j*NZ+k].tt/(u0*u0*mom[j*NZ+k].n) - SQU(mom[j*NZ+k].ut) )*fnorm;
      }
      else {
	mom[j*NZ+k].tz = mom[j*NZ+k].tr = mom[j*NZ+k].tt = 0.;
      }
      
    }
  }

  file1 = fopen(dat1, "w");
  file2 = fopen(dat2, "w");
  file3 = fopen(dat3, "w");


  for (j=0; j<nr; j++) {
    for (k=0; k<nz; k++) {  
      fprintf(file1, "% .5e % .5e % .5e \n", dr*j, dz*k, mom[j*NZ+k].tz); 
      fprintf(file2, "% .5e % .5e % .5e \n", dr*j, dz*k, mom[j*NZ+k].tr);  
      fprintf(file3, "% .5e % .5e % .5e \n", dr*j, dz*k, mom[j*NZ+k].tt); 
    }
  }
  
  fclose(file1);
  fclose(file2);
  fclose(file3);

}




void out_coords_2D( Particle pa[], size_t np, int fnorm, 
		    double omega_pe, double dz, const char *dat ) {
  FILE *file;
  double vnorm = dz/omega_pe/fnorm;
  
  file = fopen(dat, "w");
  
  for (size_t i=0; i<np; i++)  {
    fprintf(file, "% .5e % .5e % .5e % .5e % .5e \n", pa[i].p.z*dz, pa[i].p.r*dz, 
	    pa[i].p.vz*vnorm, pa[i].p.vr*vnorm, pa[i].p.vt*vnorm); 
  }
  
  fclose(file);
  
}

/***************************
 Diagnostics 
***************************/

void diagn_stability_2D( const double dens[], Particle pa[], double diag_Te[], double diag_ve[], double diag_ne[], 
			 double sign, int np, double u0, int nr, int nz, int NR, int NZ, double Omega_pe,  
			 double dr, double dz, int steps, int& check, const char *dat_err ) {

  double ne_max = 0., ne_max_old = 0.;
  double ne_av = 0., ne_av_old = 0.; // average over all k's
  int ne_max_j(-1), ne_max_k(-1), ne_av_j(-1);

  // ne/Te ratios
  double noT_max = 0., noT_max_old = 0.;
  double noT_av = 0., noT_av_old = 0.; // average over all k's
  int noT_max_j(-1), noT_max_k(-1), noT_av_j(-1);

  double fne;
  fne = SQU(1./Omega_pe);  

  double vz, vr, vt, vz2, vr2, vt2;

  FILE *file;

  // Check ne, ne/Te 
  for (int j=0; j<NR; j++) {
    for (int k=0; k<NZ; k++) {
      ne_av_old += sign*dens[j*NZ+k];
      ne_max_old = sign*dens[j*NZ+k];
      
      if ( ne_max < ne_max_old ) {
	ne_max = ne_max_old;
	ne_max_j = j;
	ne_max_k = k;
      }
    }
    if ( ne_av < ne_av_old ) {
      ne_av = ne_av_old;
      ne_av_j = j;
    }
    ne_av_old = 0.;
  }
  ne_av  *= fne/NZ; // in units of n_ref
  ne_max *= fne;    // in units of n_ref
  
  for (int j=0; j<NR; j++) {
    for (int k=0; k<NZ; k++) {
      diag_Te[j*NZ+k]=0.;
      diag_ve[j*NZ+k]=0.;
      diag_ne[j*NZ+k]=0.;
    }
  }

  for (int n=0; n<np; n++) {
    int j = (int)pa[n].p.r;
    int k = (int)pa[n].p.z;
    
    vz = pa[n].p.vz;
    vr = pa[n].p.vr;
    vt = pa[n].p.vt;
    diag_ve[j*NZ+k] += vz + vr + vt;
    
    vz2 = SQU(vz);
    vr2 = SQU(vr);
    vt2 = SQU(vt);
    diag_Te[j*NZ+k] += vz2 + vr2 + vt2;
    
    diag_ne[j*NZ+k]++;
  }
  
  for (int j=0; j<nr; j++) {
    for (int k=0; k<nz; k++) {
      diag_Te[j*NZ+k] /= diag_ne[j*NZ+k];
      //diag_Te[j*NZ+k] -= SQU(diag_ve[j*NZ+k]/diag_ne[j*NZ+k]); // don't extract main velocity
      diag_Te[j*NZ+k] /= SQU(u0);
      
      diag_ne[j*NZ+k] = sign*dens[j*NZ+k];
      
      noT_av_old += diag_ne[j*NZ+k]/diag_Te[j*NZ+k];
      noT_max_old = diag_ne[j*NZ+k]/diag_Te[j*NZ+k];
      
      if ( noT_max < noT_max_old ) {
	noT_max = noT_max_old;
	noT_max_j = j;
	noT_max_k = k;
      }
    }
    if ( noT_av < noT_av_old ) {
      noT_av = noT_av_old;
      noT_av_j = j;
    }
    noT_av_old = 0.;
  }
  noT_av  *= fne/nz; // in units of n_ref/T_ref
  noT_max *= fne;    // in units of n_ref/T_ref
  
  // Write to error file if needed
  file = fopen(dat_err, "a");
  // If error => put check == 1 => stop.
  if ( ne_av > 5. ) {
    fprintf(file, "***ERROR*** GLOBAL(%d steps): <ne>/n_ref= %.5e along j= %d \n", steps, ne_av, ne_av_j); 
    check = 1;
  }
  if ( noT_av > 2. ) {
    fprintf(file, "***ERROR*** GLOBAL(%d steps): <ne/Te>/(n_ref/T_ref)= %.5e along j= %d \n", steps, noT_av, noT_av_j);      
    check = 1;
  }
  
  // WARNINGS
  if ( ne_av > 1. ) {
    fprintf(file, "***WARNING*** GLOBAL(%d steps): <ne>/n_ref= %.5e along j= %d \n", steps, ne_av, ne_av_j); 
  }
  if ( ne_max > 5. ) {
    fprintf(file, "***WARNING*** LOCAL(%d steps): [ne]max/n_ref= %.5e at j= %d k= %d \n", steps, ne_max, ne_max_j, ne_max_k); 
  }
  if ( noT_av > 1. ) {
    fprintf(file, "***WARNING*** GLOBAL(%d steps): <ne/Te>/(n_ref/T_ref)= %.5e along j= %d \n", steps, noT_av, noT_av_j); 
  }
  if ( noT_max > 2. ) {
    fprintf(file, "***WARNING*** LOCAL(%d steps): [ne/Te]max/(n_ref/T_ref)= %.5e at j= %d k= %d \n", steps, noT_max, noT_max_j, noT_max_k); 
  }
  
  fclose(file);
}

void diagn_av_stability( const double dens[], double diag_Te[], double diag_ne[], int n_av,
			 double sign, double u0, int nr, int nz, int NR, int NZ, double Omega_pe,  
			 double dr, double dz, int steps, int *check, const char *dat_err ) {

  double ne_max = 0., ne_max_old = 0.;
  double ne_av = 0., ne_av_old = 0.; // average over all k's
  int ne_max_j(-1), ne_max_k(-1), ne_av_j(-1);

  // ne/Te ratios
  double noT_max = 0., noT_max_old = 0.;
  double noT_av = 0., noT_av_old = 0.; // average over all k's
  int noT_max_j(-1), noT_max_k(-1), noT_av_j(-1);

  if (n_av < 0) n_av = 1;
  double fne = SQU(1./Omega_pe)/n_av;  


  FILE *file;

  // Check ne, ne/Te 
  for (int j=0; j<NR; j++) {
    for (int k=0; k<NZ; k++) {
      ne_av_old += sign*dens[j*NZ+k];
      ne_max_old = sign*dens[j*NZ+k];
      
      if ( ne_max < ne_max_old ) {
	ne_max = ne_max_old;
	ne_max_j = j;
	ne_max_k = k;
      }
    }
    if ( ne_av < ne_av_old ) {
      ne_av = ne_av_old;
      ne_av_j = j;
    }
    ne_av_old = 0.;
  }
  ne_av *= fne/NZ; // in units of n_ref
  ne_max *= fne; // in units of n_ref
  
  for (int j=0; j<nr; j++) {
    for (int k=0; k<nz; k++) {
      diag_Te[j*NZ+k] /= diag_ne[j*NZ+k]; // don't extract main velocity
      diag_Te[j*NZ+k] /= SQU(u0);
      
      if ( diag_Te[j*NZ+k] > 1e-20 ) {
	noT_av_old += sign*dens[j*NZ+k]/diag_Te[j*NZ+k];
	noT_max_old = sign*dens[j*NZ+k]/diag_Te[j*NZ+k];
      }
      
      if ( noT_max < noT_max_old ) {
	noT_max = noT_max_old;
	noT_max_j = j;
	noT_max_k = k;
      }
    }
    if ( noT_av < noT_av_old ) {
      noT_av = noT_av_old;
      noT_av_j = j;
    }
    noT_av_old = 0.;
  }
  noT_av  *= fne/nz; // in units of n_ref/T_ref
  noT_max *= fne;    // in units of n_ref/T_ref

  // Write to error file if needed
  file = fopen(dat_err, "a");

  // IF ERROR => put check == 1 => stop.
  if ( ne_av > 5. ) {
    fprintf(file, "***ERROR*** GLOBAL(%d steps): <ne>/n_ref= %.5e along j= %d \n", steps, ne_av, ne_av_j); 
    *check = 1;
  }
  if ( noT_av > 2. ) {
    fprintf(file, "***ERROR*** GLOBAL(%d steps): <ne/Te>/(n_ref/T_ref)= %.5e along j= %d \n", steps, noT_av, noT_av_j);      
    *check = 1;
  }

  // WARNINGS
  if ( ne_av > 1. ) {
    fprintf(file, "***WARNING*** GLOBAL(%d steps): <ne>/n_ref= %.5e along j= %d \n", steps, ne_av, ne_av_j); 
  }
  if ( ne_max > 5. ) {
    fprintf(file, "***WARNING*** LOCAL(%d steps): [ne]max/n_ref= %.5e at j= %d k= %d \n", steps, ne_max, ne_max_j, ne_max_k); 
  }
  if ( noT_av > 1. ) {
    fprintf(file, "***WARNING*** GLOBAL(%d steps): <ne/Te>/(n_ref/T_ref)= %.5e along j= %d \n", steps, noT_av, noT_av_j); 
  }
  if ( noT_max > 2. ) {
    fprintf(file, "***WARNING*** LOCAL(%d steps): [ne/Te]max/(n_ref/T_ref)= %.5e at j= %d k= %d \n", steps, noT_max, noT_max_j, noT_max_k); 
  }
  
  fclose(file);
  
}
 

