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

  pic.h:
  Definition of structres

***********************************************************************/

#ifndef PIC_H
#define PIC_H

typedef struct  { double z, r, t; }  Vec3d;

typedef struct  { double z, r, vz, vr, vt; int m; }  PhasespaceVar;

typedef struct 	{PhasespaceVar p; }  Particle;

typedef struct  { double r; int Y; }  Sput;

//typedef struct  { double  uz, ur, ut, tz, tr, tt, eflz, eflr, 
//eflt, tzr, trt, ttz; int n; }  Moments;

typedef struct  { double  uz, ur, ut, tz, tr, tt; int n; }  Moments;


typedef struct {
  int N;
  double Emin;
  double Emax;
  double Estep;
  double Eth;
  double Wfc;
  double CS[1024];
}  Reaction;



#endif
