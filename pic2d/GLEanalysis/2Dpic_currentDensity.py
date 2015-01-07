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
# 2Dpic_currentDensity.py:
# Plots current density on anode and cathode using the 
# jhist_anode.dat and jhist_cathode.dat files
#

import sys, os, shutil
import re

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import Normalize, LogNorm
rcParams.update({'text.usetex': True})

import re

if len(sys.argv) != 5 and len(sys.argv) != 6:
    print "Usage: ./2Dpic_currentDensity.py <mintime> <maxtime> <every nth frame to analyze> <speed> [jMin:jMax]"
    exit(1)

mintime   = int(sys.argv[1])
maxtime   = int(sys.argv[2])
skipFrame = int(sys.argv[3])

speed = float(sys.argv[4])

if len(sys.argv) == 6:
    minmax = sys.argv[5].split(":")
    if len(minmax) != 2:
        print "didn't understand jMin:jMax, got '" + sys.argv[5] + "'"
    jMin = float(minmax[0])
    jMax = float(minmax[1])
else:
    jMin = -1e10
    jMax = 1e10
    print "Setting default jMin=" + str(jMin) + " [A/cm^2], jMax=" + str(jMax) + "jMax"

#Get the scaling setup
MODLOAD_parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,MODLOAD_parentdir) 
from calcScaling import InputFile
inputfile = InputFile("../input.txt")
inputfile.calcBasicParams()

#Create output folder
ofoldername = "pngs/currentDensity"
if os.path.exists(ofoldername):
    if os.path.isdir(ofoldername):
        shutil.rmtree(ofoldername)
    else:
        print "Path '" + ofoldername + "' exists, but is not a directory. Aborting!"
        exit(1)
os.mkdir(ofoldername)
print "Created directory '%s'" % (ofoldername,)

#Open input files, skip boring data
anodeFile   = open("../out/jhist_anode.dat", 'r')
assert anodeFile.readline()[:2] == "##"
match = re.match("## Anode, nr=(\d+), rmax=(\d+) ##", anodeFile.readline())
nr = match.group(1)
assert anodeFile.readline()[:2] == "##"
cathodeFile = open("../out/jhist_cathode.dat", 'r')
assert cathodeFile.readline()[:2] == "##"
assert re.match("## Cathode, nr=(\d+), rmax=(\d+) ##", cathodeFile.readline()).group(1) == nr
assert cathodeFile.readline()[:2] == "##"

nr = int(nr)
assert nr == inputfile.nr

binCenters = (np.linspace(0,nr-1,nr)+0.5)*inputfile.dZ*1e-2 #[m]
print "binCenters[um] =", binCenters*1e6
binEdges = (np.linspace(0,nr,nr+1))*inputfile.dZ*1e-2 #[m]
print "binEdges[um] =", binEdges*1e6
binEdgesSquare = []
binEdgesSquare.append(binEdges[0]);
for i in xrange(1,len(binEdges)-1):
    binEdgesSquare.append(binEdges[i]);
    binEdgesSquare.append(binEdges[i]);
binEdgesSquare.append(binEdges[-1]);
binEdgesSquare = np.asarray(binEdgesSquare)
print "binEdgesSquare[um] =", binEdgesSquare*1e6
binArea = np.pi*(binEdges[1:]**2 - binEdges[:-1]**2) #[m^2]
print "binArea[um^2] =", binArea*1e12

#Arrays for total current plot
Itime = []
Icat  = []
Iano  = []

#Arrays for current density 2D plot
j2D_time = [0.0]
j2D_cat  = []
j2D_ano  = []

outIdx = 0;
pCyclic = 0
t=0.0
tPrev = 0.0
while True:
    anoLine = anodeFile.readline()
    catLine = cathodeFile.readline()
    #print "aL = '" + anLine + "'"
    #print "cL = '" + catLine + "'"
    
    if anoLine == "":
        if catLine == "":
            print "End of files!"
            break
        else:
            print "Error, unsynchronized end of files?!?!?"
            exit(1)
    else:
        if catLine == "":
            print "Error, unsynchronized end of files?!?!?"
            exit(1)

    
    anoLineSplit = anoLine.split()
    catLineSplit = catLine.split()
    
    ts   = int(anoLineSplit[0])
    assert int(catLineSplit[0]) == ts
    t = ts*inputfile.dT*1e9 #[ns]
    integrateTime = t-tPrev
    print "integrateTime =", integrateTime, "[ns]"
    
    tPrev = t
    
    if ts < mintime:
        continue
    elif ts > maxtime:
        print "ts =", ts, "reached maxtime =", maxtime
        break
    elif pCyclic % skipFrame != 0:
        print "skipFrame", ts
        pCyclic += 1
        continue
    print "Plotting ts=", ts, ",", t, "[ns]"

    ## Plotting ##
    
    # Read the data
    catData = np.asarray(map(float, catLineSplit[1:]))
    anoData = np.asarray(map(float, anoLineSplit[1:]))
    assert len(catData) == len(anoData) == nr
    onesData = np.ones_like(catData)
    #scale to particles/bin/s
    catData  *= inputfile.N_sp/(integrateTime*1e-9)
    anoData  *= inputfile.N_sp/(integrateTime*1e-9)
    onesData *= inputfile.N_sp/(integrateTime*1e-9)
    #scale to Coulumb/m^2/s = A/m^2
    catData  *= 1.60217646e-19 / binArea
    anoData  *= 1.60217646e-19 / binArea
    onesData *= 1.60217646e-19 / binArea
    #Scale to A/cm^2
    catData  /= 1e4
    anoData  /= 1e4
    onesData /= 1e4
    #Invert old ArcPic version output format so that cathode also (mostly) positive)
    #catData *= -1

    # Plot
    def squareBins(data):
        #Convert num particles per bin into something more "plottable"
        # ( plt.step() might also work )
        dout = []
        for d in data:
            dout.append(d)
            dout.append(d)
        return np.asarray(dout);
    catDataSquare = squareBins(catData)
    anoDataSquare = squareBins(anoData)
    onesData = squareBins(onesData)
    
    plt.fill_between(binEdgesSquare*1e6, -onesData, onesData, edgecolor='gray',facecolor='gray', alpha=0.5, hatch='x')
    plt.plot(binEdgesSquare*1e6, catDataSquare, label="cathode")
    plt.plot(binEdgesSquare*1e6, anoDataSquare, label="anode")
    
    plt.ylim(jMin,jMax)
    plt.xlim(0.0,inputfile.R*1e4);
    plt.xlabel("R [$\mathrm{\mu m}$]")
    plt.ylabel("Current density [$\mathrm{A/cm^2}$]")
    plt.title("time = " + ("%.4f" % (t,) ) + " [ns]")
    plt.legend()
    #plt.savefig(ofoldername + "/currentDensity_%08d.png" %(outIdx,),dpi=300)
    
    #minCurrentDensityPossible = 1.60217646e-19*inputfile.N_sp/(integrateTime*1e-9)/np.max(binArea)
    minCurrentDensityPossible = onesData.min()
    #print "minCurrentDensityPossible", minCurrentDensityPossible, np.log10(minCurrentDensityPossible), "A/cm^2"
    #print catData
    
    plt.yscale('symlog', linthreshy=10**(int(np.log10(minCurrentDensityPossible))))
    plt.savefig(ofoldername + "/currentDensity_log_%08d.png" %(outIdx,),dpi=300)
    plt.clf()
    
    print "max current density =", max(np.max(catData),np.max(anoData)), "[A/cm^2]"
    print "min current density =", min(np.min(catData),np.min(anoData)), "[A/cm^2]"
    
    #Calculate average total current I = \int_0^nr j*(2*pi*r)*dr
    Itime.append(t)
    Icat.append( np.sum(catData*binArea*1e4) )
    Iano.append( np.sum(catData*binArea*1e4) )
    
    #Fill 2D plot
    j2D_time.append(t)
    j2D_cat.append(catData)
    j2D_ano.append(anoData)
    
    pCyclic += 1
    outIdx  += 1


#Make movie
out_dt = t/float(outIdx) #[ns]
fps = speed/out_dt
movieFileName = ofoldername + "/currentDensity.mp4"
print "Assembling movie '" + movieFileName + "' at", fps, "fps:"
print "command = '" + "ffmpeg -sameq -r "+str(fps)+" -i "+ ofoldername +"/currentDensity_%08d.png " + movieFileName + "'"
os.system("rm " + movieFileName)
os.system("ffmpeg -sameq -r "+str(fps)+" -i "+ ofoldername +"/currentDensity_%08d.png " + movieFileName)
movieFileName = ofoldername + "/currentDensity_log.mp4"
print "Assembling movie '" + movieFileName + "' at", fps, "fps:"
print "command = '" + "ffmpeg -sameq -r "+str(fps)+" -i "+ ofoldername +"/currentDensity_log_%08d.png " + movieFileName + "'"
os.system("rm " + movieFileName)
os.system("ffmpeg -sameq -r "+str(fps)+" -i "+ ofoldername +"/currentDensity_log_%08d.png " + movieFileName)


plt.figure(1)
plt.plot(Itime, Icat, label='Cathode')
plt.plot(Itime, Iano, label='Anode')
plt.xlabel("Time [ns]")
plt.ylabel("Current [A]")
plt.legend()


(j2D_time_grid,binEdges_grid) = np.meshgrid(j2D_time, binEdges*1e6)
j2D_cat = np.asarray(j2D_cat).transpose()
j2D_ano = np.asarray(j2D_ano).transpose()

plt.figure(2)
if j2D_cat.min() < 0.0:
    PC1 = plt.pcolor(j2D_time_grid,binEdges_grid, j2D_cat, norm=LogNorm(), cmap='Reds')
    CB1 = plt.colorbar(PC1)
    CB1.set_label("+j [A/cm$^2$]")

    PC2 = plt.pcolor(j2D_time_grid,binEdges_grid, -j2D_cat, norm=LogNorm(), cmap='Blues')
    CB2 = plt.colorbar(PC2)
    CB2.set_label("-j [A/cm$^2$]")

else:
    PC1 = plt.pcolor(j2D_time_grid,binEdges_grid, j2D_cat, norm=LogNorm())
    CB1 = plt.colorbar(PC1)
    CB1.set_label("-j [A/cm$^2$]")
    
plt.xlim(min(j2D_time),max(j2D_time))
plt.ylim(0.0,inputfile.R*1e4);
plt.xlabel("Time [ns]")
plt.ylabel("R [$\mathrm{\mu m}$]")
plt.title("Cathode current density [A/cm^2]")

plt.figure(3)
if j2D_ano.min() < 0.0:
    PC1 = plt.pcolor(j2D_time_grid,binEdges_grid, j2D_ano, norm=LogNorm(), cmap='Reds')
    CB1 = plt.colorbar(PC1)
    CB1.set_label("+j [A/cm$^2$]")

    PC2 = plt.pcolor(j2D_time_grid,binEdges_grid, -j2D_ano, norm=LogNorm(), cmap='Blues')
    CB2 = plt.colorbar(PC2)
    CB2.set_label("-j [A/cm$^2$]")

else:
    PC1 = plt.pcolor(j2D_time_grid,binEdges_grid, j2D_ano, norm=LogNorm())
    CB1 = plt.colorbar(PC1)
    CB1.set_label("+j [A/cm$^2$]")

plt.xlim(min(j2D_time),max(j2D_time))
plt.ylim(0.0,inputfile.R*1e4);
plt.xlabel("Time [ns]")
plt.ylabel("R [$\mathrm{\mu m}$]")
plt.title("Anode current density [A/cm^2]")

plt.show()
