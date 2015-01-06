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

  phi.cpp:
  Poisson solver using SuperLU
  Boundary conditions controlled by global variable 'BC'
  BC == 0: Phi=0 at r=nr aka infinity
  BC == 1: DPhi/Dr=0 in r=nr aka infty
  BC == 2: Phi const at r-boundaries, DPhi/Dz=0 at electrodes
  BC == 3: Periodic B.C.
  BC == 4: DPhi/Dr=0 in r=nr aka infty, alternative implementation of BC=1

  Probably the implementation of Neuman boundaries on 0-2 and 3 is flawed,
  please see the manual (CLIC note 1032)

***********************************************************************/

#include  <slu_ddefs.h>

#include  "dim.h"
#define XTRN extern
#include  "var.h"
#undef XTRN
 

void potential_factorise_2D( int nr, int nz, int NR, int NZ, double dr, double dz,
			     SuperMatrix* L, SuperMatrix* U, int** perm_c, int** perm_r ) {

  int ni, nnz;
  int *colind, *rowptr; 
  double *val;

  int *rowind, *colptr;
  double *valt;
  SuperMatrix S;
  SuperMatrix SC;
  int *etree;

  superlu_options_t options;
  SuperLUStat_t stat;
  int info;

  ni = (nr+1)*(nz+1);
  //nnz = Number of NonZeros in the matrix
  if ( BC == 0 ) {
    nnz = 2*(nr+1) + (nz-1) + 4*(nz-1) + 5*(nr-1)*(nz-1);
    //int nnz = 5*(nr+1)*(nz+1); //maximum that can be
  }
  else if ( BC == 1 ) {
    nnz = 2*(nr+1) + 3*(nz-1) + 4*(nz-1) + 5*(nr-1)*(nz-1);
  }
  else if ( BC == 4 ) {
    nnz = 2*(nr+1) + 8*(nz-1) + 5*(nr-1)*(nz-1);
  }
  else {
    printf("In potential_factorize_2D: Error with solver B.C.'s \n");
    exit(1);
  }

  printf ("In potential_factorise_2D, ni= %d, nnz= %d \n", ni, nnz);



  // I. Fill LHS
  if ( !(val = doubleMalloc(nnz)) ) ABORT("Malloc fails for val[].");
  if ( !(colind = intMalloc(nnz)) ) ABORT("Malloc fails for colind[].");
  if ( !(rowptr = intMalloc(ni+1)) ) ABORT("Malloc fails for rowptr[].");
  //fill_LHS( nr, nz, NR, NZ, dr, dz, ni, nnz, &colind, &rowptr, &val );
  //printf ("Done fill LHS. \n");
  //printf ("ni= %d, nnz= %d \n", ni, nnz);

  // LHS
  int ind=0;
  double a,b,c,e,w;
  for (int j=0; j<NR; j++) {
    for (int k=0; k<NZ; k++) {
      int i=j*NZ+k;
     
      // Cathode
      if ((i % NZ) == 0) {
	val[ind]=1.; 
	colind[ind]=i;
	rowptr[i]=ind;
	
	ind=ind+1;
      }
      
      // Anode
      else if ((i % NZ) == nz) {
	val[ind]=1.;
	colind[ind]=i;
	rowptr[i]=ind;
	
	ind=ind+1;
      }
      
      // R = infinity, vacuum
      else if (( i > nr*NZ ) && ( i < (nr*NZ+nz) )) {
	if ( BC == 0 ) {
	  val[ind]=1.;
	  colind[ind]=i;
	  rowptr[i]=ind;
	  
	  ind=ind+1;
	}
	else if ( BC == 1 ) {
	  /*
	    a=-1.; 
	    e=1.; 
	    val[ind]=a;
	    val[ind+1]=e;
	    colind[ind]=(j-1)*NZ+k;
	    colind[ind+1]=j*NZ+k;
	    rowptr[i]=ind;
	    ind=ind+2;
	  */
	  
	  // 3-p.-fla 27.8.2010
	  a=1.;
	  e=-0.75;
	  w=-0.25;
	  
	  val[ind]=w;
	  val[ind+1]=a;
	  val[ind+2]=e;
	  colind[ind]=(j-2)*NZ+k;
	  colind[ind+1]=(j-1)*NZ+k;
	  colind[ind+2]=j*NZ+k;
	  rowptr[i]=ind;
	  
	  ind=ind+3;
	}
	else if ( BC == 4) {
	  a = 2.;
	  c = 1.;
	  e = -4;
	  
	  val[ind]=a;
	  val[ind+1]=c;
	  val[ind+2]=e;
	  val[ind+3]=c;
	  colind[ind]=(j-1)*NZ+k;
	  colind[ind+1]=j*NZ+k-1;
	  colind[ind+2]=j*NZ+k;
	  colind[ind+3]=j*NZ+k+1;
	  rowptr[i]=ind;
	  
	  ind=ind+4;
	}
      }
      
      // Symmetry axis
      else if (( i > 0 ) && ( i < nz )) {
	if ( BC == 0 || BC == 1 ) {
	  b=4.; 
	  c=1.; 
	  e=-6.; 

	  val[ind]=c;
	  val[ind+1]=e;
	  val[ind+2]=c;
	  val[ind+3]=b;
	  colind[ind]=j*NZ+k-1;
	  colind[ind+1]=j*NZ+k;
	  colind[ind+2]=j*NZ+k+1;
	  colind[ind+3]=(j+1)*NZ+k;
	  rowptr[i]=ind;
	  
	  ind=ind+4;
	}
	else if (BC == 4) {
	  b=2.; 
	  c=1.; 
	  e=-4.; 

	  val[ind]=c;
	  val[ind+1]=e;
	  val[ind+2]=c;
	  val[ind+3]=b;
	  colind[ind]=j*NZ+k-1;
	  colind[ind+1]=j*NZ+k;
	  colind[ind+2]=j*NZ+k+1;
	  colind[ind+3]=(j+1)*NZ+k;
	  rowptr[i]=ind;
	  
	  ind=ind+4;	  
	}

      }
      
      // Inner parameter space
      else {
	a=1.-1./(2.*j); 
	b=1.+1./(2.*j); 
	c=1.; 
	e=-4.; 
	
	val[ind]=a;
	val[ind+1]=c;
	val[ind+2]=e;
	val[ind+3]=c;
	val[ind+4]=b;
	colind[ind]=(j-1)*NZ+k;
	colind[ind+1]=j*NZ+k-1; 
	colind[ind+2]=j*NZ+k;
	colind[ind+3]=j*NZ+k+1;
	colind[ind+4]=(j+1)*NZ+k;
	rowptr[i]=ind;
	
	ind=ind+5;
      }
    }
  }

  //Termination of rowptr
  rowptr[ni]=nnz;	 
  if (ind != nnz) {
    printf("Warning in potential_factorise_2D(): First unused ind=%d != nnz=%d", ind, nnz);
  }
  
  // II. Triangular factorisation
  //factorise( ni, nnz, &colind, &rowptr, &val, L, U, &perm_c, &perm_r );
  
  // II.1. Row -> coloumn compressed
  if ( !(rowind = intMalloc(nnz))  ) ABORT("Malloc fails for rowind[].");
  if ( !(colptr = intMalloc(ni+1)) ) ABORT("Malloc fails for colptr[].");
  if ( !(valt = doubleMalloc(nnz)) ) ABORT("Malloc fails for valt[].");
  dCompRow_to_CompCol(ni, ni, nnz, val, colind, rowptr, &valt, &rowind, &colptr);
  
  // II.2. Construct coloumn-compressed Supermatrix S
  dCreate_CompCol_Matrix(&S, ni, ni, nnz, valt, rowind, colptr, SLU_NC, SLU_D, SLU_GE);
  //dPrint_CompCol_Matrix("S",&S); // To test
  
  // II.3. Get perm_c
  set_default_options(&options);
  StatInit(&stat);
  int permc_spec = options.ColPerm;
  get_perm_c(permc_spec, &S, *perm_c);
  //print_int_vec("\nperm_c",ni,*perm_c); // To test

  // II.4. Triangular solve
  if ( !(etree = intMalloc(ni)) ) ABORT("Malloc fails for etree[].");
  sp_preorder(&options, &S, *perm_c, etree, &SC);
  int panel_size = sp_ienv(1);
  int relax = sp_ienv(2);
  dgstrf(&options, &SC, relax, panel_size, etree, NULL, 0, *perm_c, *perm_r, L, U, &stat, &info);
  //printf("Done dgstrf \n"); // To test
  //dPrint_SuperNode_Matrix("L",L);
  //dPrint_CompCol_Matrix("U",U);

  SUPERLU_FREE (val);
  SUPERLU_FREE (colind);
  SUPERLU_FREE (rowptr);
}

void potential_factorise_BC23( int nr, int nz, int NR, int NZ, double dr, double dz,
			     SuperMatrix* L, SuperMatrix* U, int** perm_c, int** perm_r ) {

  int ni, nnz;
  int *colind, *rowptr; 
  double *val;

  int *rowind, *colptr;
  double *valt;
  SuperMatrix S;
  SuperMatrix SC;
  int *etree;

  superlu_options_t options;
  SuperLUStat_t stat;
  int info;

  ni = (nr+1)*(nz+1);
  nnz = 0; //Avoid uninitialized warning
  if ( BC == 2 ) {
    nnz = 3*(nr-1) + 3*(nr-1) + (nz+1) + (nz+1) + 5*(nr-1)*(nz-1);
  }
  else if ( BC == 3 ) {
    //nnz = 4*(nr-1) + 2*(nr-1) + 2*(nz+1) + 4*(nz+1) + 5*(nr-1)*(nz-1);
    nnz = 5*(nr+1)*(nz+1); //maximum that can be
  }
  else {
    printf("Error with solver B.C.'s \n");
    exit(1);
  }
  printf ("In potential_factorise_BC23, ni= %d, nnz= %d \n", ni, nnz);
  
  // I. Fill LHS
  if ( !(val = doubleMalloc(nnz))  ) ABORT("Malloc fails for val[].");
  if ( !(colind = intMalloc(nnz))  ) ABORT("Malloc fails for colind[].");
  if ( !(rowptr = intMalloc(ni+1)) ) ABORT("Malloc fails for rowptr[].");
  //fill_LHS( nr, nz, NR, NZ, dr, dz, ni, nnz, &colind, &rowptr, &val );
  //printf ("Done fill LHS. \n");
  //printf ("ni= %d, nnz= %d \n", ni, nnz);

  // LHS
  int ind=0;
  double a,b,c,d,e,w;
  for (int j=0; j<NR; j++) {
    for (int k=0; k<NZ; k++) {
      int i=j*NZ+k;
      
      // Cathode
      if ( (k == 0) && (j > 0) && (j < nr) ) {
	if ( BC == 2 ) {
	  /*
	    d=1.;
	    e=-1.; 
	    
	    val[ind]=e;
	    val[ind+1]=d;
	    colind[ind]=j*NZ+k;
	    colind[ind+1]=j*NZ+k+1;
	    rowptr[i]=ind;
	    
	    ind=ind+2;
	  */
	  
	  // 3-p.-fla 27.8.2010
	  d=-1.;
	  e=0.75;
	  w=0.25; 
	  
	  val[ind]=e;
	  val[ind+1]=d;
	  val[ind+2]=w;
	  colind[ind]=j*NZ+k;
	  colind[ind+1]=j*NZ+k+1;
	  colind[ind+2]=j*NZ+k+2;
	  rowptr[i]=ind;
	  
	  ind=ind+3;
	}
	else if ( BC == 3 ) {
	  a=1.-1./(2.*j); 
	  b=1.+1./(2.*j); 
	  c=1.; 
	  e=-4.;
	  
	  val[ind]=a;
	  val[ind+1]=c;
	  val[ind+2]=e;
	  val[ind+3]=c;
	  val[ind+4]=b;
	  colind[ind]=(j-1)*NZ+k;
	  colind[ind+1]=j*NZ+nz-1; 
	  colind[ind+2]=j*NZ+k;
	  colind[ind+3]=j*NZ+k+1;
	  colind[ind+4]=(j+1)*NZ+k;
	  rowptr[i]=ind;
	  
	  ind=ind+5;
	}
      }
      
      // Anode
      else if ( (k == nz) && (j > 0) && (j < nr) ) {
	if ( BC == 2 ) {
	  /*
	    c=-1.; 
	    e=1.; 
	    
	    val[ind]=c;
	    val[ind+1]=e;
	    colind[ind]=j*NZ+k-1;
	    colind[ind+1]=j*NZ+k;
	    rowptr[i]=ind;
	    
	    ind=ind+2;
	  */
	  
	  // 3-p.-fla 27.8.2010
	  c=1.; 
	  e=-0.75; 
	  w=-0.25;
	  
	  val[ind]=w;
	  val[ind+1]=c;
	  val[ind+2]=e;
	  colind[ind]=j*NZ+(k-2);
	  colind[ind+1]=j*NZ+(k-1);
	  colind[ind+2]=j*NZ+k;
	  rowptr[i]=ind;
	  
	  ind=ind+3;
	}
	else if ( BC == 3 ) {
	  a=1.-1./(2.*j); 
	  b=1.+1./(2.*j); 
	  c=1.; 
	  e=-4.;
	  
	  val[ind]=a;
	  val[ind+1]=c;
	  val[ind+2]=e;
	  val[ind+3]=c;
	  val[ind+4]=b;
	  colind[ind]=(j-1)*NZ+k;
	  colind[ind+1]=j*NZ+k-1; 
	  colind[ind+2]=j*NZ+k;
	  colind[ind+3]=j*NZ+1;
	  colind[ind+4]=(j+1)*NZ+k;
	  rowptr[i]=ind;
	  
	  ind=ind+5;
	}
      }
      
      // R = infinity, vacuum
      else if ( j == nr ) {
	if ( BC == 2 ) {
	  val[ind]=1.;
	  colind[ind]=i;
	  rowptr[i]=ind;
	  
	  ind=ind+1;
	}
	else if ( BC == 3 ) {
	  /*
	    a=-1.; 
	    e=1.;
	    
	    val[ind]=a;
	    val[ind+1]=e;
	    colind[ind]=(j-1)*NZ+k;
	    colind[ind+1]=j*NZ+k;
	    rowptr[i]=ind;
	    
	    ind=ind+2;
	  */
	  
	  // 3-p.-fla 27.8.2010
	  a=1.; 
	  e=-0.75;
	  w=-0.25;
	  
	  val[ind]=w;
	  val[ind+1]=a;
	  val[ind+2]=e;
	  colind[ind]=(j-2)*NZ+k;
	  colind[ind+1]=(j-1)*NZ+k;
	  colind[ind+2]=j*NZ+k;
	  rowptr[i]=ind;
	  
	  ind=ind+3;
	}
      }
      
      // Symmetry axis
      else if ( j == 0 ) {
	if ( BC == 2 ) {
	  val[ind]=1.; 
	  colind[ind]=i;
	  rowptr[i]=ind;
	  
	  ind=ind+1;
	}
	else if ( BC == 3 ) {
	  a=4.; 
	  c=1.; 
	  e=-6.;
	  
	  val[ind]=c;
	  val[ind+1]=e;
	  val[ind+2]=c;
	  val[ind+3]=a;
	  
	  if ( k == 0 )
	    colind[ind]=j*NZ+nz-1;
	  else
	    colind[ind]=j*NZ+k-1;
	  
	  colind[ind+1]=j*NZ+k;
	  
	  if ( k == nz )
	    colind[ind+2]=j*NZ+1;
	  else
	    colind[ind+2]=j*NZ+k+1;
	  
	  colind[ind+3]=(j+1)*NZ+k;
	  rowptr[i]=ind;
	  
	  ind=ind+4;
	}
      }
      
      // Inner parameter space
      else {
	a=1.-1./(2.*j); 
	b=1.+1./(2.*j);
	c=1.; 
	e=-4.;
	
	val[ind]=a;
	val[ind+1]=c;
	val[ind+2]=e;
	val[ind+3]=c;
	val[ind+4]=b;
	colind[ind]=(j-1)*NZ+k;
	colind[ind+1]=j*NZ+k-1; 
	colind[ind+2]=j*NZ+k;
	colind[ind+3]=j*NZ+k+1;
	colind[ind+4]=(j+1)*NZ+k;
	rowptr[i]=ind;
	
	ind=ind+5;
      }
      rowptr[ni]=nnz;	 
    }
  }
  
  // II. Triangular factorisation
  //factorise( ni, nnz, &colind, &rowptr, &val, L, U, &perm_c, &perm_r );
  // II.1. Row -> coloumn compressed
  if ( !(rowind = intMalloc(nnz)) ) ABORT("Malloc fails for rowind[].");
  if ( !(colptr = intMalloc(ni+1)) ) ABORT("Malloc fails for colptr[].");
  if ( !(valt = doubleMalloc(nnz)) ) ABORT("Malloc fails for valt[].");
  dCompRow_to_CompCol(ni, ni, nnz, val, colind, rowptr, &valt, &rowind, &colptr);

  // II.2. Construct coloumn-compressed Supermatrix S
  dCreate_CompCol_Matrix(&S, ni, ni, nnz, valt, rowind, colptr, SLU_NC, SLU_D, SLU_GE);
  //dPrint_CompCol_Matrix("S",&S); // To test
  
  // II.3. Get perm_c
  set_default_options(&options);
  StatInit(&stat);
  int permc_spec = options.ColPerm;
  get_perm_c(permc_spec, &S, *perm_c);
  //print_int_vec("\nperm_c",ni,*perm_c); // To test
  
  // II.4. Triangular solve
  if ( !(etree = intMalloc(ni)) ) ABORT("Malloc fails for etree[].");
  sp_preorder(&options, &S, *perm_c, etree, &SC);
  int panel_size = sp_ienv(1);
  int relax = sp_ienv(2);
  dgstrf(&options, &SC, relax, panel_size, etree, NULL, 0, *perm_c, *perm_r, L, U, &stat, &info);
  //printf("Done dgstrf \n"); // To test
  //dPrint_SuperNode_Matrix("L",L);
  //dPrint_CompCol_Matrix("U",U);

  SUPERLU_FREE (val);
  SUPERLU_FREE (colind);
  SUPERLU_FREE (rowptr);
}

void potential_backsolve_2D( int nr, int nz, int NR, int NZ, double dz,
			     double const phi0, double const phiNz, double Phi[],
			     SuperMatrix L, SuperMatrix U, int* perm_c, int* perm_r,
			     double n_e[], double n_i[], double** rhs ) {
  
  int ni = (nr+1)*(nz+1);
  int nrhs=1;
  SuperMatrix B;

  //double phirmin, phirmax;

  SuperLUStat_t stat;
  int info;

  if (not (BC == 0 or BC == 1 or BC == 4) ) {
    printf("In potential_factorize_2D: Error with solver B.C.'s \n");
    exit(1);
  }

  // RHS vector
  for (int j=0; j<NR; j++) {
    for (int k=0; k<NZ; k++) {
      int i=j*NZ+k;
      
      // Cathode         
      if ((i % NZ) == 0) {
	(*rhs)[i]=phi0; 
      }
      
      // Anode
      else if ((i % NZ) == nz) {
	(*rhs)[i]=phiNz; 
      }
      
      // R = infinity, vacuum
      else if (( i > nr*NZ ) && ( i < (nr*NZ+nz) )) {
	if ( BC == 0 || BC == 1 ) {
	  (*rhs)[i]=0.;
	}
	else if ( BC == 4 ) {
	  (*rhs)[i]=-n_e[i]-n_i[i];
	}
      }
      
      // Symmetry axis
      else if (( i > 0 ) && ( i < nz )) {
	(*rhs)[i]=-n_e[i]-n_i[i]; // n_e<0 //CORR 29.03.2010 + -> -
      }
      
      // Inner parameter space
      else {
	(*rhs)[i]=-n_e[i]-n_i[i]; // n_e<0 //CORR 29.03.2010 + -> -
      }
    }
  }
  
  // II. Do backsolve
  //backsolve( nr, nz, Phi, rhs, L, U, perm_c, perm_r );

  // II.1. Construct RHS
  dCreate_Dense_Matrix(&B, ni, nrhs, *rhs, ni, SLU_DN, SLU_D, SLU_GE);
  //dPrint_Dense_Matrix("B",&B); // To test

  // II.2. Triangular solve
  //dPrint_SuperNode_Matrix("L_backsolve",&L); // To test
  //dPrint_CompCol_Matrix("U_backsolve",&U); // To test

  StatInit(&stat);
  trans_t trans=NOTRANS;
  dgstrs (trans, &L, &U, perm_c, perm_r, &B, &stat, &info);
  //dPrint_Dense_Matrix("B",&B); // To test

  // II.3. Return solution
  if ( info == 0 ) {
    double *sol = (double*) ((DNformat*) B.Store)->nzval;
    for (int i=0; i<ni; i++) {
      Phi[i]=sol[i];
      //printf("Phi[i] = %.5f for i = %d \n",Phi[i],i); // To test
    }
  } 
  else {
    printf("Error in SLU solver, INFO= %d\n", info);
    exit(1);
  }
  
}


void potential_backsolve_BC23( int nr, int nz, int NR, int NZ, double dz,
			     double const phi0, double const phiNz, double Phi[],
			     SuperMatrix L, SuperMatrix U, int* perm_c, int* perm_r,
			     double n_e[], double n_i[], double** rhs ) {
  
  int ni = (nr+1)*(nz+1);
  int nrhs=1;
  SuperMatrix B;
  
  double phirmin, phirmax;
  
  SuperLUStat_t stat;
  int info;

  if ( BC == 2 ) {
    phirmin = 1.;
    phirmax = 10.;
  }
  else if (BC != 3) {
    printf("In potential_backsolve_BC23: Error with solver B.C.'s \n");
    exit(1);
  }

  // I. Fill RHS
  //fill_RHS( nr, nz, NR, NZ, Phi0, PhiNz, n_e, n_i, rhs );
 
  // RHS
  int ind=0;
  for (int j=0; j<NR; j++) {
    for (int k=0; k<NZ; k++) {
      int i=j*NZ+k;
      
      // Cathode         
      if ( (k == 0) && (j > 0) && (j < nr) ) {
	if ( BC == 2 ) {
	  (*rhs)[i]=0.;
	  ind=ind+1;
	}
	else if ( BC == 3 ) {
	  (*rhs)[i]=-n_e[j*NZ]-n_i[j*NZ]-n_e[j*NZ+nz]-n_i[j*NZ+nz]; // per. density 
	  ind=ind+1;
	}
      }
      
      // Anode
      else if ( (k == nz) && (j > 0) && (j < nr) ) {
	if ( BC == 2 ) {
	  (*rhs)[i]=0.;
	  ind=ind+1;
	}
	else if ( BC == 3 ) {
	  (*rhs)[i]=-n_e[j*NZ]-n_i[j*NZ]-n_e[j*NZ+nz]-n_i[j*NZ+nz]; // per. density
	  ind=ind+1;
	}
      }
      
      // R = infinity, vacuum
      else if ( j == nr ) {
	if ( BC == 2 ) {
	  (*rhs)[i]=phirmax; 
	  ind=ind+1;
	}
	else if ( BC == 3 ) {
	  (*rhs)[i]=0.;
	  ind=ind+1;
	}
      }
      
      // Symmetry axis
      else if ( j == 0 ) {
	if ( BC == 2 ) {
	  (*rhs)[i]=phirmin; 
	  ind=ind+1;
	}
	else if ( BC == 3 ) {
	  if ( k == 0 || k == nz )
	    (*rhs)[i]=-n_e[0]-n_i[0]-n_e[nz]-n_i[nz]; // per. density
	  else
	    (*rhs)[i]=-n_e[i]-n_i[i]; 
	  ind=ind+1;
	}
      }
      
      // Inner parameter space
      else {
	(*rhs)[i]=-n_e[i]-n_i[i]; // n_e<0 //CORR 29.03.2010 + -> -
	ind=ind+1;
      }
    }
  }
  
  // II. Do backsolve
  //backsolve( nr, nz, Phi, rhs, L, U, perm_c, perm_r );

  // II.1. Construct RHS
  dCreate_Dense_Matrix(&B, ni, nrhs, *rhs, ni, SLU_DN, SLU_D, SLU_GE);
  //dPrint_Dense_Matrix("B",&B); // To test

  // II.2. Triangular solve
  //dPrint_SuperNode_Matrix("L_backsolve",&L); // To test
  //dPrint_CompCol_Matrix("U_backsolve",&U); // To test

  StatInit(&stat);
  trans_t trans=NOTRANS;
  dgstrs (trans, &L, &U, perm_c, perm_r, &B, &stat, &info);
  //dPrint_Dense_Matrix("B",&B); // To test

  // II.3. Return solution
  if ( info == 0 ) {
    double *sol = (double*) ((DNformat*) B.Store)->nzval;
    for (int i=0; i<ni; i++) {
      Phi[i]=sol[i];
      //printf("Phi[i] = %.5f for i = %d \n",Phi[i],i); // To test
    }
    
  } 
  else {
    printf("Error in SLU solver, INFO= %d\n", info);
    exit(1);
  }

}
