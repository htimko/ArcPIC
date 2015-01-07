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
# 2Dpic_density.py:
# Calculates the particle density from the r{e|Cu|Cup}...dat files
# using a histogram, and creates a density map.
# Uses equal-volume binning (not equal in DeltaR)
#

import sys, os, shutil
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm 
from matplotlib.ticker import LogFormatter

from matplotlib import gridspec

from matplotlib import rcParams
rcParams.update({'text.usetex': True})


#Get input
if len(sys.argv) != 8 and len(sys.argv) != 10:
    print "Usage: ./2Dpic_density.py <mintime> <maxtime> <every nth frame to analyse> <nBins_r> <nBins_z> <{e|Cu|Cup|qdens(L|LS|R)|all}> <speed> [<maxR[grids]> <max[grids]>]"
    exit(1)

mintime   = int(sys.argv[1])
maxtime   = int(sys.argv[2])
skipFrame = int(sys.argv[3])

nBins_r = int(sys.argv[4])
nBins_z = int(sys.argv[5])

species = sys.argv[6]
if not (species == "e" or species == "Cu" or species == "Cup" or species == "qdens" or species == "qdensL" or species == "qdensLS" or species == "qdensR" or species == "all"):
    print "species must be e, Cu, Cup, or qdensL/qdensLS/qdensR"
    exit(1)

speed = float(sys.argv[7])


#Get the scaling setup
MODLOAD_parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,MODLOAD_parentdir) 
from calcScaling import InputFile
inputfile = InputFile("../input.txt")
inputfile.calcBasicParams()

nr      = inputfile.nr  #system size [grids]
nz      = inputfile.nz  #system size [grids]

N_sp  = inputfile.N_sp  #Particle/superparticle ratio
dZ    = inputfile.dZ*1e4    #um/grid
dz    = inputfile.dz    #grid/Ldb (typically 0.5)

if len(sys.argv) == 10:
    maxR = int(sys.argv[8])
    maxZ = int(sys.argv[9])

    assert maxR <= nr
    assert maxZ <= nz
    
    nr = maxR
    nz = maxZ
    print "!!! RESTRICTING nr, nz !!!"

#Create output folder
ofoldername = "pngs/density_%s_%ix%i_%ix%i" % (species, nBins_z, nBins_r,nz,nr)
if os.path.exists(ofoldername):
    if os.path.isdir(ofoldername):
        shutil.rmtree(ofoldername)
    else:
        print "Path '" + ofoldername + "' exists, but is not a directory. Aborting!"
        exit(1)
os.mkdir(ofoldername)
print "Created directory '%s'" % (ofoldername,)

print "Got options:"
print " - mintime   =", mintime
print " - maxtime   =", maxtime
print " - skipFrame =", skipFrame
print
print " - nBins_r =", nBins_r
print " - nBins_z =", nBins_z
print " - nr      =", nr
print " - nz      =", nz
print
print " - N_sp    =", N_sp
print " - dZ      =", dZ, "[um/grid]"
print " - dz      =", dz, "[grid/Ldb]"
print
print " - species =", species
print " - speed   =", speed, "[ns/s]"
print

#Generate bin edges [dz]
print 
rEdge = np.sqrt(np.arange(nBins_r,dtype=float)+1)*(nr*dz)/np.sqrt(nBins_r)
print "Radial bin edges:"
print rEdge
#Ditto in z:
dz_bin = nz*dz/float(nBins_z) #[Ldb]
#Cell volume
Vcell = np.pi*rEdge[0]**2*dz_bin*(dZ/dz*1e-4)**3 #[cm^3]
print "Cell volume =",Vcell, "cm^3"
minDens = 1/Vcell
print "minDens =", minDens, "cm^-3"

#Generate meshgrid for plotting (bin edges)
Rgrid_points = np.asarray([0.0] + list(rEdge)) * dZ/dz
Zgrid_points = np.linspace(0,nz,nBins_z+1) * dZ
(Zgrid,Rgrid) = np.meshgrid(Zgrid_points, Rgrid_points)
#Combined plot
RgridBIG_points = np.asarray(list(-rEdge[::-1]) + [0.0] + list(rEdge)) * dZ/dz
#print RgridBIG_points
(ZgridBIG,RgridBIG) = np.meshgrid(Zgrid_points,RgridBIG_points)


def makeHisto(inFileName):

    histo = np.zeros((nBins_r,nBins_z),dtype=np.float64)
    nDropped = 0

    #Fill histogram
    inFile = open(inFileName, 'r');
    for line in inFile:
        l = line.split()
        
        z  = float(l[0]); r  = float(l[1]); #[Ldb]

        zi = int(z/dz_bin);

        ri = None
        for ri_search in xrange(len(rEdge)):
            if r <= rEdge[ri_search]:
                ri = ri_search
                break
            
        #print ri,zi, r, rEdge[ri]
        
        if zi >= nBins_z:
            #print "Warning: Dropped electron at (z,r) = (%g,%g)" %(z,r)
            nDropped += 1
            continue
        if ri == None:
            nDropped += 1
            continue
        histo[ri, zi] += 1;
    inFile.close()
    
    if nDropped > 0:
        print "Dropped", nDropped, "particles"

    #Normalize to num. physical particles / cm^3
    histo *= N_sp/Vcell;

    return histo

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

pCyclic = 0; #Used for skipping
outIdx = 0;
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
    print "Plotting ts=", ts, ",", timestamp[i], "[ns]"

    if (species == "e" or species == "Cu" or species == "Cup"):

        h = makeHisto("../out/r"+ species +stepNum[i]+".dat")
        print "histo done."
    
        #print h
        #print "log10(max):", np.log10(max(h.flatten()))

        # Plot using bin edges:
        #ticks = np.logspace(10,22,7)
        logminDens = np.floor(np.log10(1/Vcell))
        logmaxDens = max(logminDens+10,20)
        #ticks = np.power(10,np.linspace(logminDens,20,20-logminDens+1))
        ticks = np.power(10,np.linspace(logminDens,logmaxDens,logmaxDens-logminDens+1))
        if i == 0:
            print "logminDens =", logminDens, "logmaxDens =", logmaxDens, "ticks =", ticks
        plt.pcolormesh(Zgrid,Rgrid,h, norm=LogNorm(vmin=min(ticks),vmax=max(ticks)))
        plt.axis([0,max(Zgrid_points),0,max(Rgrid_points)])

        cbar = plt.colorbar(format=LogFormatter())
        cbar.set_ticks(ticks)
        cbar.set_ticklabels(map(lambda t: "$10^{"+( "%d"%(t,) )+"}$", np.log10(ticks)))
        #print ticks

        plt.xlabel("z [um]")
        plt.ylabel("r [um]")
        plt.title("Species: " + species + " , time = " + ("%.4f" % (timestamp[i],) ) + " [ns]")

    elif species.startswith("qdens"):
        hEle = makeHisto("../out/re"   +stepNum[i]+".dat")
        hCup = makeHisto("../out/rCup" +stepNum[i]+".dat")

        qDens = (hCup-hEle) #[e-/cm^3]

        if species == "qdens":
            qDens *= 1.60217646e-19 #[C/cm^3]
            #print max(qDens.flatten()), min(qDens.flatten())
            plt.pcolormesh(Zgrid,Rgrid,qDens, vmin=-5e-4,vmax=5e-4)

            cbar = plt.colorbar()
            plt.title("Charge density [$\mathrm{C}/\mathrm{cm}^3$], time = " + ("%.4f" % (timestamp[i],) ) + " [ns]")

        elif species == "qdensL":
            qDens *= 1.60217646e-19 #[C/cm^3]
            ticks = np.logspace(-8,0,9)
            plt.pcolormesh(Zgrid,Rgrid,np.abs(qDens), norm=LogNorm(vmin=min(ticks),vmax=max(ticks)))

            cbar = plt.colorbar(format=LogFormatter())
            cbar.set_ticks(ticks)
            cbar.set_ticklabels(map(lambda t: "$10^{"+( "%d"%(t,) )+"}$", np.log10(ticks)))
            plt.title("Charge density [$\mathrm{C}/\mathrm{cm}^3$], time = " + ("%.4f" % (timestamp[i],) ) + " [ns]")

        elif species == "qdensLS":
            #ticks = np.linspace(-20,20,9)
            logminDens = np.floor(np.log10(1/Vcell))
            logmaxDens = 20
            ticks = np.linspace(logminDens,logmaxDens,logmaxDens-logminDens+1)-logminDens
            ticks = np.asarray(list(-ticks[::-1]) + [0.0] + list(ticks))

            qDensL = np.log10(np.fabs(qDens))
            #print max(qDensL.flatten()), min(qDensL.flatten())

            qDensL -= logminDens-1
            for ri in xrange(nBins_r):
                for zi in xrange(nBins_z):
                    if (hEle[ri,zi] > 0 or hCup[ri,zi] > 0) and (qDensL[ri,zi] < 1 or np.isnan(qDens[ri,zi])):
                        print "zero'd ", qDens[ri,zi]
                        qDensL[ri,zi] = 0.0
            qDensL *= np.sign(qDens)
            
            # plt.get_cmap().set_bad(color='white')
            #print plt.get_cmap().get_bad()

            plt.pcolormesh(Zgrid,Rgrid,np.ma.masked_invalid(qDensL,copy=False), vmin=min(ticks),vmax=max(ticks))

            cbar = plt.colorbar()
            cbar.set_ticks(ticks)
            ticklabels = map(lambda t: "$10^{"+( "%d"%(abs(t)+logminDens,) )+"}$", ticks)
            for j in xrange(len(ticklabels)):
                if ticks[j] == 0:
                    ticklabels[j] = "$0.0$"
                elif ticks[j] < 0:
                    ticklabels[j] = "-" + ticklabels[j]
                elif ticks[j] > 0:
                    ticklabels[j] = "~" + ticklabels[j]
                
            cbar.set_ticklabels(ticklabels)
            
            plt.title("Charge density [$\mathrm{e}/\mathrm{cm}^3$], time = " + ("%.4f" % (timestamp[i],) ) + " [ns]")

        elif species == "qdensR":
            print "qdensR"
            qTot = hCup+hEle #[e-/cm^3]
            qRel = qDens/qTot
            
            plt.pcolormesh(Zgrid,Rgrid, qRel, vmin=-1.0,vmax=1.0)
            plt.colorbar()

            plt.title("Relative charge density, time = " + ("%.4f" % (timestamp[i],) ) + " [ns]")

        plt.axis([0,max(Zgrid_points),0,max(Rgrid_points)])
        plt.xlabel("z [um]")
        plt.ylabel("r [um]")
        
    elif species == "all":
        #fig = plt.figure()
        plt.title("Charge density [$\mathrm{C}/\mathrm{cm}^3$], time = " + ("%.4f" % (timestamp[i],) ) + " [ns]")

        gspec = gridspec.GridSpec(2,2)

        #left-hand plot: Charged species
        ax1 = plt.subplot(gspec[:,0])
        
        hEle = makeHisto("../out/re" +stepNum[i]+".dat")
        hCup = makeHisto("../out/rCup" +stepNum[i]+".dat")

        hCombined = np.empty((2*nBins_r,nBins_z),dtype=np.float64)
        hCombined[0:nBins_r,:] = hEle[::-1,:]
        hCombined[nBins_r:2*nBins_r,:] = hCup[:,:]

        # Plot using bin edges:
        ticks = np.logspace(10,20,6)
        plt.pcolormesh(ZgridBIG,RgridBIG,hCombined, norm=LogNorm(vmin=min(ticks),vmax=max(ticks)))
        plt.axis([0,max(Zgrid_points),-max(Rgrid_points),max(Rgrid_points)])


        plt.xlabel("z [um]")
        plt.ylabel("r [um]")
        plt.title("Cup")

        #upper-right hand plot
        ax2 = plt.subplot(gspec[0,1])

        hCu = makeHisto("../out/rCu" +stepNum[i]+".dat")
    
        # Plot using bin edges:
        ticks = np.logspace(10,20,6)
        plt.pcolormesh(Zgrid,Rgrid,hCu, norm=LogNorm(vmin=min(ticks),vmax=max(ticks)))
        plt.axis([0,max(Zgrid_points),0,max(Rgrid_points)])

        plt.title("Cu")
        plt.xlabel("z [um]")
        plt.ylabel("r [um]")

        cbar = plt.colorbar(format=LogFormatter())
        cbar.set_ticks(ticks)
        cbar.set_ticklabels(map(lambda t: "$10^{"+( "%d"%(t,) )+"}$", np.log10(ticks)))

        #Lower-right hand plot
        ax3 = plt.subplot(gspec[1,1])
        qDens = (hCup-hEle) #[e/cm^3]

        ticks = np.linspace(-8,8,9)
        plt.pcolormesh(Zgrid,Rgrid,np.log10(np.fabs(qDens))*np.sign(qDens), vmin=min(ticks),vmax=max(ticks))
        plt.axis([0,max(Zgrid_points),0,max(Rgrid_points)])
        
        cbar = plt.colorbar()
        cbar.set_ticks(ticks)
        ticklabels = map(lambda t: "$10^{"+( "%d"%(abs(t),) )+"}$", ticks)
        for i in xrange(len(ticklabels)):
            if ticks[i] < 0:
                ticklabels[i] = "-" + ticklabels[i]
        cbar.set_ticklabels(ticklabels)
        
        plt.title("qdens")
        plt.xlabel("z [um]")
        plt.ylabel("r [um]")
        
    plt.savefig(ofoldername + "/dens_%08d.png" %(outIdx,),dpi=300)
    plt.clf()

    pCyclic += 1
    outIdx += 1

movieFileName = ofoldername + "/" + os.path.basename(ofoldername) + ".mp4"
if mintime>0:
    out_dt = (timestamp[2]-timestamp[1])*skipFrame
else:
    out_dt = (timestamp[1]-timestamp[0])*skipFrame

fps = speed/out_dt
print "Assembling movie '" + movieFileName + "' at", fps, "fps:"

os.system("rm " + movieFileName)
os.system("ffmpeg -sameq -r "+str(fps)+" -i "+ ofoldername +"/dens_%08d.png " + movieFileName)
#os.system("rm _tmp*.png")
