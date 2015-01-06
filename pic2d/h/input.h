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

  input.h:
  Header file for input.cpp

***********************************************************************/

#ifndef INPUT_H
#define INPUT_H

#include <cstdio>
#include <vector>

//Toplevel input.txt parser
void input( void );

/*This reads a section on the form

(something) : <name_to_use>
*** <name1>
<Option 1> : <value> <comment bla bla>
<Option 2> : <value> <comment bla bla>
<Option 3> : <value> <comment bla bla>
///
*** <name2>
<Option 1> : <value> <comment bla bla>
(etc)
///
*** END

Arguments:
 - in_file: pointer to the input.txt file to read from
 - options_ret: Return-by-argument (see below)
 - acceptNone: Set to true to accept the wantName "None"
               corresponding to no input section.
               "None" is the returned as retName.

Return values:
 - char* containing the name of selected section
 - vector<char*> with the option values

It is the callers responsibility that char*'s are deleted.
*/
#define NAME_MAXLEN 64  //Maximum allowed length of Name
#define LINE_MAXLEN 300 //max length of line
char* readInputSection (FILE* in_file, std::vector<char*>& options_ret, bool acceptNone=false);

#endif
