#!/usr/bin/env python
#
# Copyright 2010-2014 CERN and Helsinki Institute of Physics.
# This software is distributed under the terms of the
# GNU General Public License version 3 (GPL Version 3),
# copied verbatim in the file LICENCE.md. In applying this
# license, CERN does not waive the privileges and immunities granted to it
# by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.
#
# Project website: http://arcpic.web.cern.ch/
# Developers: Helga Timko, Kyrre Sjobak
#
# This scripts reads input.txt, and calculates the system size etc.
# in physical units. Variables should have the same names as in var.h
#
# A possible future improvement would be to automatically generate dim.h
#

import re
reInt   = r"\d+"
reIntP  = "("+reInt+")"
reReal  = r"\d*\.\d*"
reRealP = "("+reReal+")"
reRealE = r"\d+\.\d+[eE]\d+"
reRealEP = "("+reRealE+")"
reNumChars = r"[\d\.eE+-]+"
reNumCharsP = "("+reNumChars+")"

import math as m
import numpy as np

#Usefull function
jFN2 = lambda(E) : 4.7133e9*((1e-9*E)**2)*np.exp(-62.338/(E*1e-9)) #A/cm^2 (V/m)

class InputFile:
    """Class which parses input.txt and stores the results internally"""
    def __init__(self, fname):
        print "InputFile parsing file '" + fname + "'..."
        ifile = open(fname, 'r');
        parseLines = None

        for line in ifile:        
            #print "LINE='"+line+"'" #Debug
            #Read lines for sections
            if parseLines != None:
                match = re.match("\*\*\*\s*END",line)
                if match:
                    if self.ArcBoundaryName != None and self.ArcBoundaryLines == None:
                        self.ArcBoundaryLines = parseLines
                    elif self.CircuitName != None and self.CircuitLines == None:
                        self.CircuitLines = parseLines
                    else:
                        print "Matched END at unknown section -- BUG!"
                        exit(1)
                    parseLines = None
                    print "Done parsing section"
                    continue
                parseLines.append(line)
                continue

            #Scaling parameters
            match = re.match("Reference density in 1/cm3:\s*" + reNumCharsP, line)
            if match:
                assert self.n_ref == None
                self.n_ref = float(match.group(1))
                print "n_ref =", self.n_ref, "[cm^-3]"
                continue
            match = re.match("Reference temperature in eV:\s*" + reNumCharsP, line)
            if match:
                assert self.T_ref == None
                self.T_ref = float(match.group(1))
                print "T_ref =", self.T_ref, "[eV]"
                continue
            match = re.match("Particles in a Debye cube Ndb:\s*" + reNumCharsP, line)
            if match:
                assert self.Ndb == None
                self.Ndb = float(match.group(1))
                print "Ndb =", self.Ndb, "[cm^-3]"
                continue

            match = re.match("Number of cells nr, nz:\s*"+ reIntP + " " + reIntP, line)
            if match:
                assert self.nr == self.nz == None
                self.nr = int(match.group(1))
                self.nz = int(match.group(2))
                print "nr =", self.nr
                print "nz =", self.nz
                continue
            match = re.match("Grid size dz in Debyes:\s*"+ reNumCharsP, line)
            if match:
                assert self.dz == None
                self.dz = float(match.group(1))
                print "dz =", self.dz
                continue
            match = re.match("Timestep dt in Omega_pe-s:\s*"+ reNumCharsP, line)
            if match:
                assert self.Omega_pe == None
                self.Omega_pe = float(match.group(1))
                print "Omega_pe =", self.Omega_pe
                continue

            match = re.match("Ion timestep dt_ion in dt-s:\s*" + reIntP, line)
            if match:
                assert self.dt_ion == None
                self.dt_ion = int(match.group(1))
                print "dt_ion =", self.dt_ion
                continue

            match = re.match("Simulation time in O_pe-s:\s*"+ reNumCharsP, line)
            if match:
                assert self.nstepsmax == None
                self.nstepsmax = int(match.group(1))
                print "nstepsmax =", self.nstepsmax, "[O_pe]"
                continue

            match = re.match("Injection timestep in dt, electrons:\s*"+ reNumCharsP, line)
            if match:
                assert self.e2inj_step == None
                self.e2inj_step = int(match.group(1))
                print "e2inj_step =", self.e2inj_step, "[dt]"
                continue

            match = re.match("Injection timestep in dt, neutrals:\s*"+ reNumCharsP, line)
            if match:
                assert self.n2inj_step == None
                self.n2inj_step = int(match.group(1))
                print "n2inj_step =", self.n2inj_step, "[dt]"
                continue

            match = re.match("Injection timestep in dt, ions:\s*"+ reNumCharsP, line)
            if match:
                assert self.i2inj_step == None
                self.i2inj_step = int(match.group(1))
                print "i2inj_step =", self.i2inj_step, "[dt]"
                continue


            #Particle boundary parameters section (initiate parsing)
            match = re.match("ArcBoundName:\s*(\S+)", line)
            if match:
                assert self.ArcBoundaryLines == None
                parseLines = []
                self.ArcBoundaryName = match.group(1)
                print "Starting to read ArcBoundName '" + self.ArcBoundaryName + "'..."
                continue
            #Circuit parameters section (initiate parsing)
            match = re.match("CircuitName:\s*(\S+)", line)
            if match:
                assert self.CircuitLines == None
                parseLines = []
                self.CircuitName = match.group(1)
                print "Starting to read CircuitName '" + self.CircuitName + "'..."
                continue
            
            #Misc
            match = re.match("mi_over_me:\s*" + reRealP, line)
            if match:
                assert self.mi_over_me == None
                self.mi_over_me = float(match.group(1))
                print "mi_over_me =", self.mi_over_me
                continue
        
        ifile.close()
        if parseLines != None:
            print "Ran off the end of file while reading parseLines - file syntax error or BUG"
            exit(1)
        print "Sucessfully parsed file."
    ### END def __init__()
    
    #Extra parsing functions
    def parseArcBoundaryParams(self):
        print
        print "Parsing and calculating ArcBoundaryParams:"
        relevantLines = self.getRelevantSection(self.ArcBoundaryLines,self.ArcBoundaryName)
        self.ArcBoundaryParams = {}
        for line in relevantLines:
            match = re.match("beta\s*:\s*"+ reNumCharsP, line)
            if match:
                assert not ("beta" in self.ArcBoundaryParams)
                self.ArcBoundaryParams["beta"] = float(match.group(1))
                print "beta =", self.ArcBoundaryParams["beta"]
                try:
                    self.ArcBoundaryParams["beta_E"] = self.CircuitParams["E"]*self.ArcBoundaryParams["beta"]
                    print "beta_E =", self.ArcBoundaryParams["beta_E"], "[V/cm] =", self.ArcBoundaryParams["beta_E"]*1e-7, "[GV/m]"
                    self.ArcBoundaryParams["beta_E_J"] = jFN2(self.ArcBoundaryParams["beta_E"]*1e2)
                    print "beta_E_J =", self.ArcBoundaryParams["beta_E_J"], "[A/cm^2]"
                except:
                    print "Couldn't calculate beta_E and/or beta_E_J"
                continue
            match = re.match("ring_beta\s*:\s*"+ reNumCharsP, line)
            if match:
                assert not ("ring_beta" in self.ArcBoundaryParams)
                self.ArcBoundaryParams["ring_beta"] = float(match.group(1))
                print "ring_beta =", self.ArcBoundaryParams["ring_beta"]
                try:
                    self.ArcBoundaryParams["ring_beta_E"] = self.CircuitParams["E"]*self.ArcBoundaryParams["ring_beta"]
                    print "ring_beta_E =", self.ArcBoundaryParams["ring_beta_E"], "[V/cm] =", self.ArcBoundaryParams["ring_beta_E"]*1e-7, "[GV/m]"
                    self.ArcBoundaryParams["ring_beta_E_J"] = jFN2(self.ArcBoundaryParams["ring_beta_E"]*1e2)
                    print "ring_beta_E_J =", self.ArcBoundaryParams["ring_beta_E_J"], "[A/cm^2] (NOT INCLUDING ALPHA!)"
  
                except all as e:
                    print e
                    print "Couldn't calculate ring_beta_E"
                continue
            l = line
            if l.endswith("\n"):
                l = l[:-1]
            print "Un-matched line: '" + l + "'"

    def parseCircuitParams(self):
        print
        print "Parsing and calculating CircuitParams:"
        #Parsing
        relevantLines = self.getRelevantSection(self.CircuitLines,self.CircuitName)
        self.CircuitParams = {}
        for line in relevantLines:
            match = re.match("U0\s*:\s*"+ reNumCharsP, line)
            if match:
                assert not ("U0" in self.CircuitParams)
                self.CircuitParams["U0"] = float(match.group(1))
                print "U0 =", self.CircuitParams["U0"]
                continue
            match = re.match("UNz\s*:\s*"+ reNumCharsP, line)
            if match:
                assert not ("UNz" in self.CircuitParams)
                self.CircuitParams["UNz"] = float(match.group(1))
                print "UNz =", self.CircuitParams["UNz"]
                continue
            match = re.match("UNz_ampl\s*:\s*"+ reNumCharsP, line)
            if match:
                assert not ("UNz_ampl" in self.CircuitParams)
                self.CircuitParams["UNz_ampl"] = float(match.group(1))
                print "UNz_ampl =", self.CircuitParams["UNz_ampl"]
                continue

            l = line
            if l.endswith("\n"):
                l = l[:-1]
            print "Un-matched line: '" + l + "'"
        #Calculate derived quantities
        if not "U0" in self.CircuitParams:
            self.CircuitParams["U0"] = 0.0
            print "U0 =", self.CircuitParams["U0"], "(by default value)"

        if "UNz" in self.CircuitParams and "U0" in self.CircuitParams:
            self.CircuitParams["deltaU"] = self.T_ref*(self.CircuitParams["UNz"]-self.CircuitParams["U0"])
        elif "UNz_ampl" in self.CircuitParams:
            self.CircuitParams["deltaU"] = self.T_ref*self.CircuitParams["UNz_ampl"]
        print "deltaU =", self.CircuitParams["deltaU"], "[V] (cathode rel. anode)"
        self.CircuitParams["E"] = self.CircuitParams["deltaU"]/self.Z
        print "E =", self.CircuitParams["E"], "[V/cm] =", self.CircuitParams["E"]*1e-9, "[GV/cm], =", self.CircuitParams["E"]*1e-7, "[GV/m]"

    def getRelevantSection(self,lines,secName):
        "Return a list of strings corresponding to the correctly named section from the input lines"
        ret = None
        done = False
        for line in lines:
            match = re.match("\*\*\*\s*("+secName+")\s*",line)
            if match and line[match.end(1):-1].strip() == '':
                print "Found matching section '" + secName +"' in line '" + line[:-1] + "'"
                assert ret == None
                assert done == False
                ret = []
                continue
            if ret != None and done == False:
                if line.strip() == "///":
                    print "End of the section."
                    done = True
                    continue
                ret.append(line)
        if ret == None:
            print "ERROR in getRelevantSection: Couldn't find section '" + secName + "'"
            exit(1)
        if done == False: #ret != None if we are here
            print "ERROR in getRelevantSection: Ran off the end while parsing section '" + secName + "'"
            exit(1)
        return ret
    ### END getRelevantSection()

    ## Input parameters ##

    #Scaling parameters:
    n_ref     = None #Reference density in cm^-3. In code n_e = n_ref.
    T_ref     = None #Reference temperature in eV. In code T_e = T_ref.
    Ndb       = None #Superparticles in a Debye cube at reference density

    nr        = None #Number of cells, r direction
    nz        = None #Number of cells, z direction
    dz        = None #Grid size in Debye lengths
    Omega_pe  = None #Time step in omega_pe^-1

    dt_ion    = None #Ion timestep as multiple of normal timesteps
    
    nstepsmax = None #Simulation time in Omega_pe^-1

    e2inj_step = None #How often [dt] to inject electrons
    n2inj_step = None #How often [dt] to inject neutrals
    i2inj_step = None #How often [dt] to inject ions

    #Parameters for the arc boundary
    ArcBoundaryName   = None
    ArcBoundaryLines  = None
    ArcBoundaryParams = None #Dict with parsed CircuitParameters,
                             # see function parseBoundaryParams() for info

    #Parameters for the circuit
    CircuitName   = None 
    CircuitLines  = None
    CircuitParams = None #Dict with parsed CircuitParameters,
                         # see function parseCircuitParams() for info
    
    #Misc parameters
    mi_over_me = None #H/electron mass ratio

    ## Calculated parameters ##
    Ldb   = None #Debye length in cm
    O_pe  = None #Plasma frequency in s^-1

    dZ    = None #Grid size in cm
    dT    = None #Time step in s
    Z     = None #System size in cm (Z-direction)
    R     = None #System size in cm (R-direction)

    N_sp  = None #Particle/superparticle ratio

    def calcBasicParams(self):
        """Calculates system size and scalings"""
        print
        print "Calculated basic parameters:"

        self.Ldb = 7.43e2*m.sqrt(self.T_ref/self.n_ref) #NOTE: This is the actual constant in code,
                                                        # not exactly equal sqrt(552635) from doc
        print "Ldb =", self.Ldb, "[cm] =", self.Ldb*1e7, "[nm]"
        
        self.O_pe = 56414.6*m.sqrt(self.n_ref)
        print "O_pe =", self.O_pe, "[s^-1]"

        self.dT = self.Omega_pe/self.O_pe
        print "dT =", self.dT, "[s] =", self.dT*1e15, "[fs]"

        self.dZ = self.Ldb*self.dz
        print "dZ =", self.dZ, "[cm] =", self.dZ*1e4, "[um]"

        self.R = self.dZ*self.nr
        print "R =", self.R, "[cm] =", self.R*1e4, "[um]"
        
        self.Z = self.dZ*self.nz
        print "Z =", self.Z, "[cm] =", self.Z*1e4, "[um]"
        
        self.T = self.nstepsmax/self.O_pe
        print "T =", self.T, "[s] = ", self.T*1e9, "[ns]"

        self.N_sp = self.n_ref*self.Ldb**3/self.Ndb
        print "N_sp =", self.N_sp
    ### END def calcBasicParams()

### END class InputFile

if __name__ == "__main__":
    inputfile = InputFile("input.txt");
    inputfile.calcBasicParams()
    inputfile.parseCircuitParams()
    inputfile.parseArcBoundaryParams()

    print
    print "Memory parameters for dim.h:"
    #TODO...
    # +Also print rought memory requirements
    # +Autogenerate dim.h
    print "\t Remember to set NGR = nr+1 =", inputfile.nr+1
    print "\t Remember to set NGZ = nz+1 =", inputfile.nz+1
