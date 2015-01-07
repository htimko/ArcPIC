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

  engy.cpp:
  Energy balance calucation 

***********************************************************************/

#include "pic.h"
#include "mydef.h"

#include <cstddef>


void energy_electrons_2D( Particle el[], size_t ne, double *We )
{
  double vte2;
  extern double v_te;

  
  *We = 0.;
  vte2 = 2.*v_te*v_te;

  for (size_t n=0; n<ne; n++)
    {
      *We += SQU(el[n].p.vz) + SQU(el[n].p.vr) + SQU(el[n].p.vt);
    }
  *We /= vte2;		/* We/Te, multiply by 1D: 1/2, 2D: 2/2, 3D: 3/2 */

}




void energy_all_2D( Particle el[], size_t ne, double *We, double *Wf, 
		    const double E_r[], const double E_z[], int NR, int NZ, 
		    double Ndb, double omega_pe, double dz)
{
     double vte2; 
     extern double v_te;
    

     *We = 0.;
     *Wf = 0.;

     vte2 = 2.*v_te*v_te;	

     for (size_t n=0; n<ne; n++)
     	{
	  *We += SQU(el[n].p.vz) + SQU(el[n].p.vr) + SQU(el[n].p.vt);
     	}
     *We /= vte2;	/* We/Te */    
     

     for (int j=0; j<NR; j++)  
       {
	 for (int k=0; k<NZ; k++)
	   {  
	     *Wf += SQU(E_z[j*NR+k]) + SQU(E_r[j*NR+k]);	
	   }
     }
     *Wf *= 4.*(Ndb*dz/SQU(omega_pe))/vte2;	/* Wf/Te, 4*q/Te, Var. E[] enth. E/2 */	

}




void energy_total_2D( Particle el[], size_t ne, Particle ion[], size_t ni, Particle neu[], size_t nn,
		      double *We, double *Wi, double *Wn, double *Wf, double *Wt,
		      double me_over_MCup, double me_over_MCu,
		      const double E_r[], const double E_z[], int NR, int NZ, 
		      double Vcell[], double Ndb, double omega_pe, double dz )
{
     double vte2; 
     extern double v_te;
    

     *We = 0.; // electrons
     *Wi = 0.; // ions
     //     *Wn = 0.; // neutrals
     *Wf = 0.; // E-field
     *Wt = 0.; // total energy

     vte2 = 2.*v_te*v_te;	

     for (size_t n=0; n<ne; n++)
     	{
	  *We += SQU(el[n].p.vz) + SQU(el[n].p.vr) + SQU(el[n].p.vt);
     	}
     *We /= vte2;	// We/Te     

     for (size_t n=0; n<ni; n++)
     	{
	  *Wi += SQU(ion[n].p.vz) + SQU(ion[n].p.vr) + SQU(ion[n].p.vt);
     	}
     *Wi /= me_over_MCup*vte2;	// We/Te   
     
     for (size_t n=0; n<nn; n++)
     	{
	  *Wn += SQU(neu[n].p.vz) + SQU(neu[n].p.vr) + SQU(neu[n].p.vt);
     	}
     *Wn /= me_over_MCu*vte2;	// We/Te   
     

     for (int j=0; j<NR; j++)  
       {
	 for (int k=0; k<NZ; k++)
	   {  
	     // correct with Vcell!!!
	     *Wf += (SQU(E_z[j*NR+k]) + SQU(E_r[j*NR+k])) * Vcell[j]; // only 2-dim. contribution!!!
	   }
     }
     *Wf *= 4.*(Ndb*dz*SQU(dz/omega_pe))/vte2;   // Wf/Te, 4*q/Te, Var. E[] enth. E/2 

     
     *Wt = *We + *Wi + *Wn + *Wf; 
     
}
 



void kin_pot_en( Particle el[], size_t ne, Particle ion[], size_t ni, Particle neu[], size_t nn,
		 double *We, double *Wi, double *Wn, double *Wp, double *Wt,
		 double me_over_MCup, double me_over_MCu,
		 const double Phi[], int NR, int NZ, double omega_pe, double dz )
{
     int j, k;
     double vte2;
     double hr,hz,Vr;
     extern double v_te;
    

     *We = 0.; // electrons kin en
     *Wi = 0.; // ions kin en
     *Wn = 0.; // neutrals kin en
     *Wp = 0.; // potential energy
     *Wt = 0.; // total energy

     vte2 = 2.*v_te*v_te;	

     for (size_t n=0; n<ne; n++)
     	{
	  *We += SQU(el[n].p.vz) + SQU(el[n].p.vr) + SQU(el[n].p.vt);
     	}
     *We /= vte2;	// We/Te     

     for (size_t n=0; n<ni; n++)
     	{
	  *Wi += SQU(ion[n].p.vz) + SQU(ion[n].p.vr) + SQU(ion[n].p.vt);
     	}
     *Wi /= me_over_MCup*vte2;	// Wi/Te   
     
     for (size_t n=0; n<nn; n++)
     	{
	  *Wn += SQU(neu[n].p.vz) + SQU(neu[n].p.vr) + SQU(neu[n].p.vt);
     	}
     *Wn /= me_over_MCu*vte2;	// Wn/Te   
     

     for (size_t n=0; n<ne; n++)
     	{
	  hr  = el[n].p.r;
	  j   = (int)hr;
	  hr -= j;       

	  hz  = el[n].p.z;
	  k   = (int)hz;
	  hz -= k;   

	  Vr  = (1-hr)*(1-hz)*Phi[j*NZ+k] + (1-hr)*hz*Phi[j*NZ+(k+1)] + hr*(1-hz)*Phi[(j+1)*NZ+k] + hr*hz*Phi[(j+1)*NZ+(k+1)]; 
	  
	  *Wp -= Vr;
	}

     for (size_t n=0; n<ni; n++)
     	{
	  hr  = ion[n].p.r;
	  j   = (int)hr;
	  hr -= j;       

	  hz  = ion[n].p.z;
	  k   = (int)hz;
	  hz -= k;   

	  Vr  = (1-hr)*(1-hz)*Phi[j*NZ+k] + (1-hr)*hz*Phi[j*NZ+(k+1)] + hr*(1-hz)*Phi[(j+1)*NZ+k] + hr*hz*Phi[(j+1)*NZ+(k+1)]; 
	  
	  *Wp += Vr;
	}
     *Wp *= SQU(dz/omega_pe);  // Wp/Te


    
     *Wt = *We + *Wi + *Wn + *Wp; 
     
}
 
