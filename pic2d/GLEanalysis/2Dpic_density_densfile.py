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
# 2Dpic_density_densfiles.py:
# Calculates the density using the density output files from ArcPIC2D.
# These gives the density at gridpoints using the Verbonceur volumes,
# averaged over n_ave_time steps.
#

import sys, os, shutil
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from matplotlib.ticker import LogFormatter, LogLocator, MaxNLocator
import matplotlib.colorbar as colorbar

from matplotlib import gridspec

from matplotlib import rcParams,rc
rcParams.update({'text.usetex': True}) #slow
DPI = 500

#Get input
if len(sys.argv) < 5 or len(sys.argv) > 8:
    print "Usage: ./2Dpic_density_densfiles.py <mintime> <maxtime> <every nth frame to analyse> <{e|Cu|Cup|qdens|all|all2|allRow}> {<axisShape='square'|'image'>} {<FPS=5>} {cutR=MAX|number}"
    exit(1)

mintime   = int(sys.argv[1])
maxtime   = int(sys.argv[2])
skipFrame = int(sys.argv[3])

species = sys.argv[4]
if not (species == "e" or species == "Cu" or species == "Cup" or species == "qdens" or species == "all" or species == "all2" or species == "allRow"):
    print "species must be e, Cu, Cup, qdens, or all(2|Row)"
    exit(1)

if len(sys.argv) >= 6:
    axisShape = sys.argv[5]
    if not (axisShape == "square" or axisShape=="image"):
        print "axisShape must be 'square' or 'image'"
        exit(1)
else:
    axisShape = 'square'
if len(sys.argv) >= 7:
    fps = float(sys.argv[6])
else:
    fps = 5

if len(sys.argv) >= 8:
    cutR = sys.argv[7]
    if cutR != "MAX":
        try:
            cutR = float(cutR)
        except ValueError:
            print "cutR must be 'MAX' or a number of um where the R-axis should be cut"
            exit(1)
    else:
        cutR = None
else:
    cutR == None

#Get the scaling setup
MODLOAD_parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,MODLOAD_parentdir) 
from calcScaling import InputFile
inputfile = InputFile("../input.txt")
inputfile.calcBasicParams()

#Create output folder
ofoldername = "pngs/density_densfiles_%s" % (species,)
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
print " - species   =", species
print " - axisShape =", axisShape
print " - fps       =", fps
print

#Read IV curves
if species == "all2":
    (circ_ts, circ_deltaQ, circ_U) = np.loadtxt("../circuit.dat", usecols=(0,1,2), unpack=True)

    convFactor2 = 1.60217657e-19*inputfile.N_sp/inputfile.dT # superparticles/timestep -> Amps
    convFactor3 = (inputfile.dz/inputfile.Omega_pe)**2*inputfile.T_ref #dimless potential -> Volts

    circ_t = circ_ts*inputfile.dT*1e9  #ns
    circ_I = circ_deltaQ * convFactor2 #A
    circ_U = circ_U * convFactor3      #V

def plotIV(time=None):
    ax1 = plt.gca()
    ax1.plot(circ_t, circ_I, "b-", label="Circuit current")
    ax1.set_ylabel("Current [A]", color='b')
    ax1.set_xlabel("t [ns]")
    for tl in ax1.get_yticklabels():
        tl.set_color('b')
    
    ax2 = ax1.twinx()
    ax2.plot(circ_t, circ_U, "r-", label="Circuit current")
    ax2.set_ylabel("Voltage [V]", color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    
    if time:
        ax1.axvline(x=time, color="k",ls="--")
    

#Function to read and rescale files
def readDensFile(fname):
    (rList, zList, densList) = np.loadtxt(fname,unpack=True)
    rList *= inputfile.Ldb*1e4 #um
    zList *= inputfile.Ldb*1e4 #um
    densList *= inputfile.n_ref #cm^-3
    return (np.reshape( rList,(inputfile.nr+1,-1) ),
            np.reshape( zList,(inputfile.nr+1,-1) ),
            np.reshape(densList,(inputfile.nr+1,-1) ) )

#Colorbar label?
drawCB  = True

#Function to plot a single species into an existing subplot
def plot_singleSpecies(r,z,dens):
    CS0 = plt.contourf(z,r,dens,100);
    if drawCB:
        CB0 = plt.colorbar(CS0)
        CB0.set_label("Density [cm$^{-3}$]")
    plt.xlabel("z [um]")
    plt.ylabel("r [um]")
    if axisShape == "image":
        plt.axis('image')
    if cutR != None:
        plt.ylim(0,cutR)
    return CS0
def plot_singleSpecies_log(r,z,dens):
    if dens.max() > 0:            
        CS0 = plt.contourf(z,r,dens,locator=LogLocator(), vmin=1e10, vmax=1e20)
        if drawCB:
            CB0 = plt.colorbar(CS0)
            CB0.set_label("Density [cm$^{-3}$]")
    else:
        print " ! Plotting empty ! "
        CS0 = plt.contourf(z,r,dens, vmin=1e10, vmax=1e20, colors='white')
        if drawCB:
            CB0 = plt.colorbar(CS0,ticks=[])
            CB0.set_label("Density [cm$^{-3}$]")
    plt.xlabel("z [\\textmu m]")
    plt.ylabel("r [\\textmu m]")
    if axisShape == "image":
        plt.axis('image')
    if cutR != None:
        plt.ylim(0,cutR)
    return CS0

def plot_qdens(r,z,dens):
    CS0 = plt.contourf(z, r, dens)
    if drawCB:
        CB0 = plt.colorbar(CS0)
        CB0.set_label("Charge density [cm$^{-3}$]")
    plt.xlabel("z [um]")
    plt.ylabel("r [um]")
    if axisShape == "image":
        plt.axis('image')
    return CS0
def plot_qdens_log(r,z,dens, dens0):
    CS0 = plt.contourf(z,r,dens0, colors='k') #mark areas with data but zero density
    if dens.max() > 0.0:
        CS1 = plt.contourf(z, r, dens,norm=LogNorm(vmin=1e10, vmax=1e20))
        if drawCB:
            CB1 = plt.colorbar(CS1)
    else:
        CS1 = None
        if drawCB:
            CB1 = plt.colorbar(CS0, ticks=[])
    
    if dens.min() < 0.0:
        CS2 = plt.contourf(z, r, -dens,norm=LogNorm(vmin=1e10, vmax=1e20), hatches='x')
        if drawCB:
            CB2 = plt.colorbar(CS2)
    else:
        CS2 = None
        if drawCB:
            CB2 = plt.colorbar(CS0, ticks=[])
    
    if drawCB:
        CB1.set_label("Positive charge density [cm$^{-3}$]")
        CB2.set_label("Negative charge density [cm$^{-3}$]")
    plt.xlabel("z [um]")
    plt.ylabel("r [um]")
    if axisShape == "image":
        plt.axis('image')
    if cutR != None:
        plt.ylim(0,cutR)
    return (CS0, CS1, CS2)

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

timeOutFile = open(ofoldername+"/timeOutFile.txt",'w')
timeOutFile.write("# timestep, t[ns], outIdx \n")
timeOutFile.flush()

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
    print "Plotting ts=", ts, ",", timestamp[i], "[ns], outIdx=", outIdx
    timeOutFile.write(str(ts) + " " + str(timestamp[i]) + " " + str(outIdx) + "\n")
    timeOutFile.flush()

    if (species == "e" or species == "Cu" or species == "Cup"):
        if species == "e":
            foo = "n"+species
        else:
            foo = species
        (r,z,dens) = readDensFile("../out/" + foo + stepNum[i] + ".dat")

        # plot_singleSpecies(r,z,dens)
        # plt.title("Density [cm$^{-3}$] of " + species + ", time = %.3f [ns]" % (timestamp[i],))
        # plt.savefig(ofoldername + "/dens_%08d.png" %(outIdx,),dpi=DPI)
        # plt.clf()

        plot_singleSpecies_log(r,z,dens)
        plt.title("Density [cm$^{-3}$] of " + species + ", time = %.3f [ns]" % (timestamp[i],))
        plt.savefig(ofoldername + "/dens_log_%08d.png" %(outIdx,),dpi=DPI)
        plt.clf()
    
    elif species == "qdens":
        (r_e,z_e,dens_e)       = readDensFile("../out/ne"  + stepNum[i] + ".dat")
        (r_Cup,z_Cup,dens_Cup) = readDensFile("../out/Cup" + stepNum[i] + ".dat")
        dens = dens_Cup - dens_e
        
        plot_qdens(r_e, z_e, dens)
        plt.title("Density [e/cm$^{3}$], time = %.3f [ns]" % (timestamp[i],))
        plt.savefig(ofoldername + "/dens_%08d.png" %(outIdx,),dpi=DPI)
        plt.clf()

        #True where one of dens is > 0, else true
        dens0_bool = np.logical_or(dens_Cup > 0.0, dens_e > 0.0)
        #0.0 where there are no data, masked if there are
        dens0 = np.ma.masked_array(np.zeros_like(dens), dens0_bool)        
        plot_qdens_log(r_e, z_e, dens, dens0)
        plt.title("Density [e/cm$^{3}$], time = %.3f [ns]" % (timestamp[i],))
        plt.savefig(ofoldername + "/dens_log_%08d.png" %(outIdx,),dpi=DPI)
        plt.clf()

    elif species == "all" or species == "all2" or species == "allRow":
        drawCB = False
        #load data
        (r_e,z_e,dens_e)       = readDensFile("../out/ne"  + stepNum[i] + ".dat")
        (r_Cup,z_Cup,dens_Cup) = readDensFile("../out/Cup" + stepNum[i] + ".dat")
        (r_Cu,z_Cu,dens_Cu)    = readDensFile("../out/Cu"  + stepNum[i] + ".dat")

        if species!="allRow":
            (fig, axes) = plt.subplots(nrows=2, ncols=2)
            titletext_Y = 0.92
        else:
            textwidth = 0.98*449.40363*1.0/72.27 #inches, for 2Dpic paper
            (fig, axes) = plt.subplots(nrows=1, ncols=3,figsize=(textwidth,textwidth/3.5))
            #plt.setp(axes.flat, aspect=1.0, adjustable='box-forced')
            titletext_Y = 0.90
            
            rc('font',**{'family':'serif','serif':['Times'],'size':10})
                        
        #neutrals
        if species!="allRow":
            ax1 = plt.subplot(2,2,1)
        else:
            ax1 = plt.subplot(1,3,1)
        Cu_CS0 = plot_singleSpecies_log(r_Cu,z_Cu,dens_Cu)
        #plt.title("Neutrals")
        plt.text(0.5, titletext_Y, "Neutrals",
                 horizontalalignment='center',
                 fontsize=10,
                 transform = ax1.transAxes)

        if species!="allRow":
            ax2 = plt.subplot(2,2,2)
        else:
            ax2 = plt.subplot(1,3,2,sharey=ax1)
        Cup_CS0 = plot_singleSpecies_log(r_Cup,z_Cup,dens_Cup)
        #plt.title("Ions")
        plt.text(0.5, titletext_Y, "Ions",
                 horizontalalignment='center',
                 fontsize=10,
                 transform = ax2.transAxes)

        if species!="allRow":
            ax3 = plt.subplot(2,2,3)
        else:
            ax3 = plt.subplot(1,3,3,sharey=ax1)
        e_CS0 = plot_singleSpecies_log(r_e,z_e,dens_e)
        #plt.title("Electrons")
        plt.text(0.5, titletext_Y, "Electrons",
                 horizontalalignment='center',
                 fontsize=10,
                 transform = ax3.transAxes)

        if species != "allRow":
            ax4 = plt.subplot(2,2,4)
            if species == "all":
                qdens = dens_Cup - dens_e
                dens0_bool = np.logical_or(dens_Cup > 0.0, dens_e > 0.0)
                dens0 = np.ma.masked_array(np.zeros_like(qdens), dens0_bool)
                (qdens_CS0, qdens_CS1, qdens_CS2) = plot_qdens_log(r_e,z_e,qdens,dens0)
                #plt.title("Charge density")
                plt.text(0.5, titletext_Y, "Charge",
                         horizontalalignment='center',
                         fontsize=10,
                         transform = ax4.transAxes,
                         color='red')
            elif species == "all2":
                plotIV(time=timestamp[i])
                plt.gca().xaxis.set_major_locator(MaxNLocator(nbins=5))
        else:
            plt.setp( ax2.get_yticklabels(), visible=False)
            ax2.set_ylabel("")
            #ax2.set_xlabel("z [\textmu m]")
            plt.setp( ax3.get_yticklabels(), visible=False)
            ax3.set_ylabel("")
            #ax3.set_xlabel("z [\textmu m]")
        ### COLORBARS AND ADJUSTMENTS ###
        if species == "all":
            fig.subplots_adjust(right=0.75)
            cax = fig.add_axes([0.8, 0.1, 0.03, 0.8])
        elif species == "all2":
            fig.subplots_adjust(right=0.8)
            cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        elif species == "allRow":
            fig.subplots_adjust(right=0.87, left=0.055, top=0.96, bottom=0.2, wspace=0.07)
            cax = fig.add_axes([0.885, 0.2, 0.03, 0.76])

        my_levels = np.logspace(10,20,11)
        #print my_levels
        CBpos = colorbar.ColorbarBase(cax,cmap=Cu_CS0.cmap,norm=Cu_CS0.norm, boundaries=my_levels, values=my_levels)

        if species == "all":
            cax2 = fig.add_axes([0.9, 0.1, 0.03, 0.8])

            if qdens_CS2 != None: #Plot empty cbar in case of no negative charge density
                CBneg = colorbar.colorbar_factory(cax2,qdens_CS2)
                CBneg.values = my_levels
                CBneg.boundaries = my_levels
                CBneg.set_ticks(np.logspace(10,20,6))
                # CBneg.update_ticks()
                CBneg.config_axis()
                CBneg.draw_all()
        if species != "allRow":
            plt.suptitle("Densities, time = %.3f [ns]" % (timestamp[i],))
        else:
            CBpos.set_label("Particle density [cm\\textsuperscript{-3}]")

        plt.savefig(ofoldername + "/dens_log_%08d.png" %(outIdx,),dpi=DPI)
        plt.savefig(ofoldername + "/dens_log_%08d.pdf" %(outIdx,))
        plt.clf()
        
    pCyclic += 1
    outIdx += 1

timeOutFile.close()

if species != "all" and species != "all2":
    movieFileName = ofoldername + "/dens.mp4"
    print "Assembling movie '" + movieFileName + "' at", fps, "fps:"
    os.system("rm " + movieFileName)
    command = "ffmpeg -sameq -r "+str(fps)+" -i "+ ofoldername +"/dens_%08d.png " + movieFileName
    print "Command ='" + command + "'"
    #os.system(command)
    print "Skipping .mp4!"
    movieFileName = ofoldername + "/dens.gif"
    command = "convert " +ofoldername +"/dens_*.png -layers Optimize " + movieFileName
    print "Command ='" + command + "'"
    #os.system(command)
    print "Skipping convert!"

movieFileName = ofoldername + "/dens_log.mp4"
print "Assembling movie '" + movieFileName + "' at", fps, "fps:"
os.system("rm " + movieFileName)
command = "ffmpeg -sameq -r "+str(fps)+" -i "+ ofoldername +"/dens_log_%08d.png " + movieFileName
print "Command = '" + command + "'"
#os.system(command)
print "Skipping .mp4!"
print "Converting to .gif:"
movieFileName = ofoldername + "/dens_log.gif"
command = "convert " +ofoldername +"/dens_log_*.png -layers Optimize -delay " + str(100/fps) + " " + movieFileName
print "Command = '" + command + "'"
#os.system(command)
print "Skipping convert!"
