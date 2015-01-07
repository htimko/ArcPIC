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

  my_time.cpp:
  Reads clock for outputting

***********************************************************************/

#include  <time.h> 
#include <stddef.h>

time_t my_time(void)
{
  time_t     stime;
  struct tm  *ltime;
   
  stime = time(NULL);
  ltime = localtime(&stime);
   
  /*ltime->tm_hour += 2; */
  return  mktime(ltime);
}




double now_msec()
{
  double time_msec;
  
  time_msec = (double)clock() / (double)CLOCKS_PER_SEC *1.e6;

  return  time_msec; 
}
