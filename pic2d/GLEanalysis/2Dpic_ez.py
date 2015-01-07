#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2010-2015 CERN and Helsinki Institute of Physics.
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
# 2Dpic_ez.py:
# Plots Ez on the cathode as a function of r
#


import sys, os, shutil
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm 
#from matplotlib.ticker import LogFormatter

#from matplotlib import gridspec

from matplotlib import rcParams
rcParams.update({'text.usetex': True})


#Get input
if len(sys.argv) != 4 and len(sys.argv) != 5:
    print "Usage: ./2Dpic_ez.py <mintime> <maxtime> <every nth frame to analyse> [<yrange>]"
    exit(1)

mintime   = int(sys.argv[1])
maxtime   = int(sys.argv[2])
skipFrame = int(sys.argv[3])

ymin = -5000
ymax = 300
if len(sys.argv) == 5:
    yrange = sys.argv[4].split(":")
    if len(yrange) != 2:
        print "yrange syntax: 'min:max'"
    ymin = float(yrange[0])
    ymax = float(yrange[1])

#Get the scaling setup
MODLOAD_parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,MODLOAD_parentdir) 
from calcScaling import InputFile
inputfile = InputFile("../input.txt")
inputfile.calcBasicParams()

rFactor     = 1e4*inputfile.Ldb #output units -> um
fieldFactor = 510.998928e3/2.99792458e10**2 * inputfile.Ldb * inputfile.O_pe**2 * 1e-4 #output units -> MV/m

#Create output folder
ofoldername = "pngs/efield_z_cathode"
if os.path.exists(ofoldername):
    if os.path.isdir(ofoldername):
        print "Removing " + ofoldername
        shutil.rmtree(ofoldername)
    else:
        print "Path '" + ofoldername + "' exists, but is not a directory. Aborting!"
        exit(1)
os.mkdir(ofoldername)
print "Created directory '%s'" % (ofoldername,)

print "Got options:"
print " - mintime     =", mintime
print " - maxtime     =", maxtime
print " - skipFrame   =", skipFrame
print
print " - yrange      = " + str(ymin) + ":" + str(ymax)
print
print " - rFactor     =", rFactor
print " - fieldFactor =", fieldFactor
print

#Get list of output files and timestamp
stepNum   = []
timestamp = [] #[ns]
timeIdx = open("../out/timeIndex.dat")
timeIdx.readline() #Skip first line
if mintime == 0:
    stepNum.append("00000000")
    timestamp.append(0.0)
for line in timeIdx:
    ls = line.split()
    stepNum.append(ls[0])
    timestamp.append(float(ls[1]))
timeIdx.close()

allFieldList = []

### Time step loop ###
pCyclic = 0; #Used for skipping
outIdx = 0;
firstTime = None
finalTime = None
for i in xrange(len(stepNum)):
    ts = int(stepNum[i])
    
    if ts < mintime:
        continue
    elif ts > maxtime:
        print "ts =", ts, "reached maxtime =", maxtime
        break
    elif pCyclic % skipFrame != 0:
        print "skipFrame", ts
        pCyclic += 1
        continue
    fname = "../out/Ez" + stepNum[i] + ".dat"
    print "Plotting ts=", ts, ",", timestamp[i], "[ns], fname='" + fname + "' "

    #Read data for this timestep
    fieldFile = open(fname, 'r')
    rList = []
    #zList = []
    fieldList = []
    for line in fieldFile:
        l = line.split()
        if len(l) == 0:
            continue
        z = float(l[1])
        if (z > 0.0):
            #Want field at cathode only
            continue
        rList.append(float(l[0]))
        #zList.append(z)
        fieldList.append(float(l[2]))
    fieldFile.close()
    rList = np.asarray(rList)*rFactor
    fieldList = np.asarray(fieldList)*fieldFactor
    allFieldList.append(fieldList);
    
    plt.plot(rList,fieldList)
    plt.ylim(ymin,ymax)

    plt.xlabel("Radius [um]")
    plt.ylabel("$\mathrm{E_z}$ [MV/m]")
    plt.title("Cathode field profile, time=%.3f ns" %(timestamp[i],))

    plt.savefig(ofoldername + "/fieldProfile_%08d.png" %(outIdx,),dpi=300)
    plt.clf()
    
    if firstTime == None:
        firstTime = timestamp[i]
    finalTime = timestamp[i]
    
    pCyclic += 1
    outIdx += 1    

plt.figure(1)
#allFieldList = np.rot90(np.asarray(allFieldList))
allFieldList = np.transpose(np.asarray(allFieldList))
plt.imshow(allFieldList, aspect='auto', origin='lower',extent=[firstTime,finalTime, 0,inputfile.R*1e4])
plt.colorbar()

plt.xlabel("Time [ns]")
plt.ylabel("Radius [um]")
plt.title("Cathode field $\mathrm{E_z}$ [MV/m]")

plt.savefig(ofoldername + "/fieldTimeMap.png")

plt.figure(2)
plt.imshow(-allFieldList, aspect='auto', origin='lower',extent=[firstTime,finalTime, 0,inputfile.R*1e4])
plt.colorbar()

plt.xlabel("Time [ns]")
plt.ylabel("Radius [um]")
plt.title("Cathode field $\mathrm{E_z}$ [MV/m]")

plt.savefig(ofoldername + "/fieldTimeMap_inv.png")

plt.show()
