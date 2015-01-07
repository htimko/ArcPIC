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

  outp.h:
  File names for output, most set by calling file_names_2D()

***********************************************************************/

#define LEN_FILENAME 28

XTRN char fn_e[LEN_FILENAME], fn_i[LEN_FILENAME], fn_n[LEN_FILENAME];
XTRN char fphi[LEN_FILENAME];
XTRN char fEz[LEN_FILENAME], fEr[LEN_FILENAME];
XTRN char fv_ez[LEN_FILENAME],   fv_er[LEN_FILENAME],   fv_et[LEN_FILENAME];
XTRN char fT_ez[LEN_FILENAME],   fT_er[LEN_FILENAME],   fT_et[LEN_FILENAME];
XTRN char fv_iz[LEN_FILENAME],   fv_ir[LEN_FILENAME],   fv_it[LEN_FILENAME];
XTRN char fT_iz[LEN_FILENAME],   fT_ir[LEN_FILENAME],   fT_it[LEN_FILENAME];
XTRN char fr_e[LEN_FILENAME],    fr_i[LEN_FILENAME],    fr_n[LEN_FILENAME]; 
XTRN char fvdf_ez[LEN_FILENAME], fvdf_er[LEN_FILENAME], fvdf_eabs[LEN_FILENAME];
XTRN char fvdf_iz[LEN_FILENAME], fvdf_ir[LEN_FILENAME], fvdf_iabs[LEN_FILENAME];
XTRN char fvdf_nz[LEN_FILENAME], fvdf_nr[LEN_FILENAME], fvdf_nabs[LEN_FILENAME];

XTRN char frestart[LEN_FILENAME], fbckp[LEN_FILENAME];

XTRN char ferr[LEN_FILENAME];

//Also have file from Arcbounds::writefile_removedParticles()
