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
# 2Dpic_potentialMap.py:
# Plots potential in gap as function of r, z at different time steps
#

import sys, os, shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.ticker import Locator, MultipleLocator
import matplotlib.colorbar

class MinorSymLogLocator(Locator):
    """
    Dynamically find minor tick positions based on the positions of
    major ticks for a symlog scaling.
    """
    def __init__(self, linthresh):
        """
        Ticks will be placed between the major ticks.
        The placement is linear for x between -linthresh and linthresh,
        otherwise its logarithmically
        """
        self.linthresh = linthresh

    def __call__(self):
        'Return the locations of the ticks'
        majorlocs = self.axis.get_majorticklocs()

        # iterate through minor locs
        minorlocs = []

        # handle the lowest part
        for i in xrange(1, len(majorlocs)):
            majorstep = majorlocs[i] - majorlocs[i-1]
            if abs(majorlocs[i-1] + majorstep/2) < self.linthresh:
                ndivs = 10
            else:
                ndivs = 9
            minorstep = majorstep / ndivs
            locs = np.arange(majorlocs[i-1], majorlocs[i], minorstep)[1:]
            minorlocs.extend(locs)

        return self.raise_if_exceeds(np.array(minorlocs))

    def tick_values(self, vmin, vmax):
        raise NotImplementedError('Cannot get tick locations for a '
                                  '%s type.' % type(self))

#from matplotlib import gridspec

from matplotlib import rcParams,rc
rcParams.update({'text.usetex': True})
DPI = 500
rcParams.update({'savefig.dpi':DPI})
rc('font',**{'family':'serif','serif':['Times'],'size':8})
textwidth = 449.40363*1.0/72.27 #inches, for 2Dpic paper
FIGSIZE  = (0.5*0.95*textwidth, 0.5*0.98*textwidth/1.61803398875) #Other plots
#FIGSIZE1 = (0.5*0.95*textwidth, 0.5*0.95*textwidth) #potential plots
FIGSIZE1 = (0.23*textwidth, 0.23*textwidth) #potential plots (half area)
FIGSIZE2 = (0.95*textwidth, 0.1*0.95*textwidth/1.61803398875) #colorbar
FIGSIZE3  = (0.98*textwidth, 0.5*0.98*textwidth/1.61803398875) #1D selected pot plot

#Get input
if len(sys.argv) != 4 and len(sys.argv) != 5 and len(sys.argv) != 6:
    print "Usage: ./2Dpic_potentialMap.py <mintime> <maxtime> <every nth frame to analyse> [potVmin:potVmax] [fps]"
    exit(1)

mintime   = int(sys.argv[1])
maxtime   = int(sys.argv[2])
skipFrame = int(sys.argv[3])

potVmin = None
potVmax = None
if len(sys.argv) == 5 or len(sys.argv) == 6:
    minmax = sys.argv[4].split(":")
    if len(minmax) != 2:
        print "didn't understand potVmin:potVmax, got '" + sys.argv[4] + "'"
        exit(1)
    potVmin = float(minmax[0])
    potVmax = float(minmax[1])

fps = 10
if len(sys.argv) == 6:
    fps = float(sys.argv[5])

#Get the scaling setup
MODLOAD_parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,MODLOAD_parentdir) 
from calcScaling import InputFile
inputfile = InputFile("../input.txt")
inputfile.calcBasicParams()

rFactor   = 1e4*inputfile.Ldb #output units -> um
potFactor = inputfile.T_ref   #output units -> Volts

#Create output folder
ofoldername = "pngs/potentialMap"
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
print " - potVmin     =", potVmin
print " - potVmax     =", potVmax
print
print " - fps         =", fps
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

potLevels = np.arange(potVmin,potVmax, 25) #25 V / division
print "potLevels=", potLevels
firstContours = None

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
    fname = "../out/phi" + stepNum[i] + ".dat"
    print "Plotting ts=", ts, ",", timestamp[i], "[ns], fname='" + fname + "' "
    
    #Read data for this timestep
    potFile = open(fname, 'r')
    rList = []
    zList = []
    potList = []
    for line in potFile:
        l = line.split()
        if len(l) == 0:
            continue
        rList.append( float(l[0]) )
        zList.append( float(l[1]) )
        potList.append(float(l[2]))
    potFile.close()
    rList = np.asarray(rList)*rFactor
    zList = np.asarray(zList)*rFactor
    potList = np.asarray(potList)*potFactor
    
    rList = np.reshape(rList,(inputfile.nr+1,-1))
    zList = np.reshape(zList,(inputfile.nr+1,-1))
    potList = np.reshape(potList,(inputfile.nr+1,-1))
    
    plt.figure(1,figsize=FIGSIZE1,dpi=DPI)
    plt.clf()
    
    contours = plt.contour(zList,rList,potList, potLevels, vmin=potVmin, vmax=potVmax,linewidths=0.5);
    plt.setp(contours.collections[4], linewidth=1.5)
    #plt.colorbar().set_label("Potential [V]")
    #plt.clabel(contours)
    if outIdx == 0:
        firstContours = contours
    plt.xlabel("z [\\textmu m]")
    plt.ylabel("r [\\textmu m]")
    # plt.title("Potential [V], time = %.3f" % (timestamp[i],))

    plt.savefig(ofoldername + "/potential_%08d.png" %(outIdx,),dpi=DPI)
    
    plt.axis('image')
    plt.ylim(0,6.0)
    #plt.subplots_adjust(right=0.97, left=0.1, top=0.99, bottom=0.1)
    plt.subplots_adjust(right=0.99, left=0.16, top=0.96, bottom=0.21)
    #plt.show()#For debug of subplots_adjust
    plt.savefig(ofoldername + "/potential_zoom_%08d.png" %(outIdx,),dpi=DPI)

    plt.figure(3,figsize=FIGSIZE,dpi=DPI)
    plt.clf()
    
    #plot potential as function of z along different r's
    for j in reversed(xrange(inputfile.nr+1)):
        plt.plot(zList[j,:],potList[j,:])

    plt.xlabel("z [um]")
    plt.ylabel("V [V]")
    plt.title("Potential [V], time = %.3f, view from axis" % (timestamp[i],))

    plt.xlim(0.0,inputfile.Z*1e4)
    plt.ylim(potVmin, potVmax)

    plt.savefig(ofoldername + "/potential1D_%08d.png" %(outIdx,),dpi=DPI)
    #plt.clf()

    #plot potential as function of z along different r, somewhat spaced and zoomed
    # for j in reversed(xrange(0,200,10)):
    #     plt.plot(zList[j,:],potList[j,:])

    # plt.xlim(0.0,1.0)
    # plt.ylim(potVmin, potVmax)

    # plt.savefig(ofoldername + "/potential1D_zoom_%08d.png" %(outIdx,),dpi=DPI)
    # plt.clf()
    
    plt.figure(4,figsize=FIGSIZE3,dpi=DPI)
    plt.clf()

    #plt.plot(zList[0],potList[0,:], label="0.0 \\textmu m")
    #plt.plot(zList[9],potList[9,:], label="0.5 \\textmu m")
    plt.plot(zList[19],potList[19,:], label="1.0 \\textmu m", ls="-")
    #plt.plot(zList[29],potList[29,:], label="1.5 \\textmu m")
    plt.plot(zList[39],potList[39,:], label="2.0 \\textmu m", ls="--")
    #plt.plot(zList[49],potList[49,:], label="2.5 \\textmu m")
    plt.plot(zList[59],potList[59,:], label="3.0 \\textmu m", ls=":")
    #plt.plot(zList[69],potList[69,:], label="3.5 \\textmu m")

    plt.legend(loc=9,frameon=False, fontsize=8, ncol=6)

    #plt.subplots_adjust(right=0.99, left=0.07, top=0.98, bottom=0.08)
    plt.subplots_adjust(right=0.99, left=0.07, top=0.96, bottom=0.17)
    plt.xlabel("z [\\textmu m]")
    plt.ylabel("Electric potential [V]")
    plt.xlim(0,6)
    plt.ylim(0,1e4)
    plt.yscale('symlog',linthreshy=10)
    #plt.minorticks_on()
    plt.gca().yaxis.set_minor_locator(MinorSymLogLocator(10))
    #plt.show()
    plt.savefig(ofoldername + "/potential1D_selected_%08d.png" %(outIdx,),dpi=DPI)

    pCyclic += 1
    outIdx += 1    

#colorbar
# plot potential contours
colorbar_fig = plt.figure(2,dpi=DPI,figsize=FIGSIZE2)
colorbar_ax = colorbar_fig.add_axes([0.025, 0.40, 0.95, 0.5])
#colorbar_cmap = matplotlib.colors.Normalize(vmin=potVmin,vmax=potVmax)
#matplotlib.colorbar.ColorbarBase(colorbar_ax,orientation='horizontal',values=potLevels,norm=colorbar_cmap)
colorbar = matplotlib.colorbar.Colorbar(colorbar_ax,firstContours,orientation='horizontal')
colorbar_linewidths = colorbar.ax.get_children()[4].get_linewidth()
colorbar_linewidths[4] = 4
colorbar.ax.get_children()[4].set_linewidth(colorbar_linewidths)
#plt.colorbar(contours).set_label("Potential [V]")
plt.savefig(ofoldername + "/potential_colorbar.png",dpi=DPI)

#plt.show()


#Movie assembly
def movieAssembly(movieName, imageBasename, fps):

    movieFileName = ofoldername + "/" + os.path.basename(movieName) + ".mp4"
    if mintime>0:
        out_dt = (timestamp[2]-timestamp[1])*skipFrame
    else:
        out_dt = (timestamp[1]-timestamp[0])*skipFrame

    print "Assembling movie '" + movieFileName + "' at", fps, "fps:"
    
    ffmpegCommand = "ffmpeg -sameq -r "+str(fps)+" -i "+ ofoldername +"/" + imageBasename + "_%08d.png " + movieFileName
    print "Command: '" + ffmpegCommand + "'"
    
    os.system("rm " + movieFileName)
    os.system(ffmpegCommand)

    #Make some space after ffmpeg output
    print
    print

#movieAssembly("potential","potential", fps)
#movieAssembly("potential1D","potential1D", fps)
