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
# 2Dpic_temperature.py:
# Calculates the particle temperature from the r{e|Cu|Cup}...dat files.
#

import sys, os, shutil
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm 
from matplotlib.ticker import LogFormatter

from matplotlib import gridspec

from matplotlib import rcParams,rc
rcParams.update({'text.usetex': True})
DPI = 500
rcParams.update({'savefig.dpi':DPI})
rc('font',**{'family':'serif','serif':['Times'],'size':8})
textwidth = 449.40363*1.0/72.27 #inches, for 2Dpic paper
FIGSIZE = (0.5*0.95*textwidth, 0.5*0.95*textwidth/1.61803398875)

#Get input
if len(sys.argv) != 11:
    print "Usage: ./2Dpic_temperature.py <mintime> <maxtime> <every nth frame to analyse> <min_r> <max_r> <min_z> <max_z> <species=e|Cu|Cup> <speed> <tempunit=K|eV>"
    exit(1)

mintime   = int(sys.argv[1])
maxtime   = int(sys.argv[2])
skipFrame = int(sys.argv[3])

min_r = float(sys.argv[4])
max_r = float(sys.argv[5])

min_z = float(sys.argv[6])
max_z = float(sys.argv[7])

species = sys.argv[8]
if not (species == "e" or species == "Cu" or species == "Cup"):
    print "species must be one of 'e', 'Cu', 'Cup'"
    exit(1)

speed = float(sys.argv[9])

tempunit = sys.argv[10]
if not (tempunit == "K" or tempunit == "eV"):
    print "tempunit must be one of 'K' or 'eV'"
    exit(1)

#Get the scaling setup
MODLOAD_parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,MODLOAD_parentdir) 
from calcScaling import InputFile
inputfile = InputFile("../input.txt")
inputfile.calcBasicParams()

nr      = inputfile.nr  #system size [grids]
nz      = inputfile.nz  #system size [grids]

N_sp  = inputfile.N_sp  #Particle/superparticle ratio
dZ    = inputfile.dZ    #um/grid
dz    = inputfile.dz    #Ldb/grid (typically 0.5)

Ldb   = inputfile.Ldb   #Debye length [cm]
O_pe  = inputfile.O_pe  #Plasma angular frequency [s^-1]

if max_r > inputfile.R*1e4:
    print
    print "************************************************************************************************"
    print "   max_r = %f [um] can't be bigger than system, setting to system size R = %f [um]" % (max_r, inputfile.R*1e4)
    print "************************************************************************************************"
    print
    max_r = inputfile.R*1e4
if max_z > inputfile.Z*1e4:
    print
    print "************************************************************************************************"
    print "   max_z = %f [um] can't be bigger than system, setting to system size Z = %f [um]" % (max_z, inputfile.Z*1e4)
    print "************************************************************************************************"
    print
    max_z = inputfile.Z*1e4

#Create output folder
ofoldername = "pngs/temperature_%s_minR=%.1f_maxR=%.1f_minZ=%.1f_maxZ=%.1f" % (species, min_r, max_r, min_z, max_z)
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
print " - min_r =", min_r, "[um]"
print " - max_r =", max_r, "[um]"
print " - min_z =", min_z, "[um]"
print " - max_z =", max_z, "[um]"
print
print " - species =", species
print
print " - speed   =", speed, "[ns/s]"
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

rFactor = dZ/dz*1e4 # output format -> um
vFactor = Ldb*O_pe * 1e-2; # output format -> m/s

plotYmax = 0; #y axis scale should grow monotonically


plotXmin = 0
plotXmax = 0
if species == "e":
    plotXmin = -10e6
    plotXmax = 30e6
elif species == "Cu":
    plotXmin = -0.05e6
    plotXmax = 0.05e6
elif species == "Cup":
    plotXmin = -0.05e6
    plotXmax =  0.05e6
plotXrange = np.linspace(plotXmin, plotXmax, 200)

circuitFile = open("../circuit.dat", 'r'); #Used to get the bias voltage
circuitFile.readline(); circuitFile.readline(); #Skip headers

#prepare for summary plot over time
log_time = []
log_Tz   = []
log_Tr   = []
log_Tt   = []
log_N    = []
log_dens = []

### Time step loop ###
pCyclic = 0; #Used for skipping
outIdx = 0;
fig = plt.figure(dpi=DPI,figsize=FIGSIZE)
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
    fname = "../out/r" + species + stepNum[i] + ".dat"
    print "Plotting ts=", ts, ",", timestamp[i], "[ns], fname='" + fname + "' ",
    
    #Read the particle file
    particleFile = open(fname,'r')
    VZ = []; muZ = 0.0; sigmaZ = 0.0;
    VR = []; muR = 0.0; sigmaR = 0.0;
    VT = []; muT = 0.0; sigmaT = 0.0;
    for line in particleFile:
        l = line.split()
        z = float(l[0])*rFactor # [um]
        r = float(l[1])*rFactor # [um]

        if not (z > min_z and z < max_z and r > min_r and r < max_r):
            continue        

        vz = float(l[2])*vFactor # [m/s]
        vr = float(l[3])*vFactor # [m/s]
        vt = float(l[4])*vFactor # [m/s]
        
        VZ.append(vz)
        VR.append(vr)
        VT.append(vt)
        
        muZ += vz;
        sigmaZ += vz*vz

        muR += vr
        sigmaR += vr*vr

        muT += vt;
        sigmaT += vt*vt

    N = len(VZ)
    print "N=", N

    if N == 0:
        #Todo: Rather plot an empty histogram?
        print "num particles = 0, SKIP!"
        continue

    #Read the circuit file
    voltage = 0;
    for line in circuitFile:
        #print ts, line[:-1]
        l = line.split()
        circTS = int(l[0])
        if circTS < ts:
            continue
        assert circTS == ts
        
        Udimless = float(l[2])
        voltage = inputfile.T_ref*Udimless*(inputfile.dz/inputfile.Omega_pe)**2
        break
    mass = 0;
    #Put theoretical max velocity on the chart
    vmax_theor = None
    if species == "e":
        voltage *= max_z/(inputfile.Z*1e4)
        mass = 9.10938291e-31 #kg
        vmax_theor = np.sqrt(2*voltage*1.60217657e-19/mass) 
        plt.axvline(vmax_theor, ls=':', color='k')
    elif species == "Cu":
        mass = (63.546*inputfile.mi_over_me/1.0073-1)*9.10938291e-31
    elif species == "Cup": 
        mass = (63.546*inputfile.mi_over_me/1.0073-1)*9.10938291e-31
        #Anyway off-the-chart
        #plt.axvline(-np.sqrt(2*voltage*1.60217657e-19/mass), ls='--', label="theo. max")

    
    #Postprocess distributions
    if N > 1:
        muZ /= N
        sigmaZ = np.sqrt( (sigmaZ - N*muZ*muZ)/(N-1) )
        tempZ = sigmaZ**2 * mass / 1.3806488e-23
        #print muZ, sigmaZ, tempZ

        muR /= N
        sigmaR = np.sqrt( (sigmaR - N*muR*muR)/(N-1) )
        tempR = sigmaR**2 * mass / 1.3806488e-23
        #print muR, sigmaR, tempR

        muT /= N
        sigmaT = np.sqrt( (sigmaT - N*muT*muT)/(N-1) )
        tempT = sigmaT**2 * mass / 1.3806488e-23    
        #print muT, sigmaT, tempT
    else:
        muZ /= N
        sigmaZ = 0.0
        tempZ = 0.0

        muR /= N
        sigmaR = 0.0
        tempR = 0.0

        muT /= N
        sigmaT = 0.0
        tempT = 0.0

    nBins = int(np.sqrt(np.ceil(N/5.0)))
    binWidth = (plotXmax-plotXmin)/float(nBins)
    print "nBins = ", nBins, "binwidth = ", binWidth*1e-6, "range = ", plotXmax*1e-6, plotXmin*1e-6, "[Mm/s]"
    norm = N*inputfile.N_sp
    histWeights = np.ones(N)*inputfile.N_sp/binWidth#*1e7;
    

    histVZ = plt.hist(np.asarray(VZ), nBins, weights=histWeights, range=(plotXmin,plotXmax), label='$v_z$', histtype="stepfilled", alpha=0.3, color='blue')
    plt.plot(plotXrange, norm * np.exp(-(plotXrange-muZ)**2/(2*sigmaZ**2))/(sigmaZ*np.sqrt(2*np.pi)),
             color='blue', ls='-.')

    histVR = plt.hist(np.asarray(VR), nBins, weights=histWeights, range=(plotXmin,plotXmax), label='$v_r$', histtype="stepfilled", alpha=0.3, color='green')
    plt.plot(plotXrange, norm * np.exp(-(plotXrange-muR)**2/(2*sigmaR**2))/(sigmaR*np.sqrt(2*np.pi)),
             color='green', ls='-.')

    histVT = plt.hist(np.asarray(VT), nBins, weights=histWeights, range=(plotXmin,plotXmax), label='$v_\\theta$', histtype="stepfilled", alpha=0.3,color='red')
    plt.plot(plotXrange, norm * np.exp(-(plotXrange-muT)**2/(2*sigmaT**2))/(sigmaT*np.sqrt(2*np.pi)),
             color='red', ls='-.')
    
    plt.gca().get_xaxis().get_major_formatter().set_powerlimits((-3, 3))
    plt.gca().get_yaxis().get_major_formatter().set_powerlimits((-3, 3))
    
    #For testing normalization
    # testData = np.random.random(N)*6e7-1e7
    # print "Average velocity density =", N*inputfile.N_sp/6e7*1e6, "particles/Mm/s"
    # plt.hist(testData, nBins, weights=histWeights, range=(plotXmin,plotXmax), histtype="stepfilled", alpha=0.3,color='orange')
    
    
    plt.xlim(plotXmin, plotXmax)

    if plt.ylim()[1] == 1.0:
        Ylim = max((max(histVZ[0]),max(histVR[0]),max(histVT[0])))
        Ylim = 10**np.ceil(np.log10(Ylim))
        plt.ylim(0,Ylim)
        print "Setting Ylim =", Ylim, "< 1.0"
    #Grow Y axis monotonically with time
    if plt.axis()[3] > plotYmax:
        plotYmax = plt.axis()[3]
    plt.ylim(0,plotYmax);

    textXoffset = 0.05
    textYoffset = 0.9
    if species == "e":
        textYoffset = 0.52
        textXoffset = 0.7

        (Ylim_min, Ylim_max) = plt.ylim()
        if vmax_theor == vmax_theor: #Protect against NaN
            plt.text(vmax_theor, 0.75*Ylim_max, r"$\sqrt{2 V_\mathrm{gap} / m_e}$",
                     horizontalalignment='right',
                     verticalalignment='center',
                     rotation='vertical')
        
        plt.text(0.9, 0.1,r"e\textsuperscript{-}",
                 horizontalalignment='center',
                 verticalalignment='center',
                 transform = plt.gca().transAxes,
                 bbox=dict(facecolor='gray', alpha=0.3))
    elif species == "Cup":
        plt.text(0.9, 0.1,r"Cu\textsuperscript{+}",
                 horizontalalignment='center',
                 verticalalignment='center',
                 transform = plt.gca().transAxes,
                 bbox=dict(facecolor='gray', alpha=0.3))
    elif species == "Cu":
        plt.text(0.9, 0.1,r"Cu",
                 horizontalalignment='center',
                 verticalalignment='center',
                 transform = plt.gca().transAxes,
                 bbox=dict(facecolor='gray', alpha=0.3))
    
    eVtempFactor = 1
    eVorKstring = "K"
    if tempunit == "eV":
        eVtempFactor = 1.3806488e-23 / 1.60217657e-19
        eVorKstring = "eV"

    # plt.text(textXoffset, textYoffset,"t = %.4f ns" % (timestamp[i]),
    #          horizontalalignment='left',
    #          verticalalignment='center',
    #          transform = plt.gca().transAxes)
    plt.text(textXoffset, textYoffset,"$T_z$ = %.3g %s" % (tempZ*eVtempFactor, eVorKstring),
             horizontalalignment='left',
             verticalalignment='center',
             transform = plt.gca().transAxes)
    plt.text(textXoffset, textYoffset-0.08,"$T_r$ = %.3g %s" % (tempR*eVtempFactor, eVorKstring),
             horizontalalignment='left',
             verticalalignment='center',
             transform = plt.gca().transAxes)
    plt.text(textXoffset, textYoffset-0.16,"$T_\\theta$ = %.3g %s" % (tempT*eVtempFactor, eVorKstring),
             horizontalalignment='left',
             verticalalignment='center',
             transform = plt.gca().transAxes)

    #print "<T> [eV] =", (tempZ+tempR+tempT)/3.0 * 1.3806488e-23 / 1.60217657e-19

    plt.text(textXoffset, textYoffset-0.24,"N = %.3g" % (N*inputfile.N_sp),
             horizontalalignment='left',
             verticalalignment='center',
             transform = plt.gca().transAxes)

    volume = (max_z-min_z)*np.pi*(max_r**2-min_r**2)*1e-12 #[cm^3]
    density = inputfile.N_sp*N/volume # [particles / cm^3]
    # plt.text(textXoffset, textYoffset-0.25,"$\\langle \\rho \\rangle$ = %.3g cm$^{-3}$" % (density),
    #          horizontalalignment='left',
    #          verticalalignment='center',
    #          transform = plt.gca().transAxes)

    log_time.append(timestamp[i])
    log_Tz.append(tempZ)
    log_Tr.append(tempR)
    log_Tt.append(tempT)
    log_N.append(N*inputfile.N_sp)
    log_dens.append(density)


    plt.legend(loc=1,ncol=1,fontsize=8,frameon=False) #rather put this in the text.
    # plt.title("Velocity distribution for " + species +
    #           ", region (Z,R) [um]: (%.1f,%.1f)-(%.1f,%.1f)" % (min_z, min_r, max_z, max_r))

    plt.xlabel("Velocity component [m/s]")
    #plt.ylabel(r"Particles / m/s ")
    plt.ylabel(r"$\mathrm{d}N/\mathrm{d}v \left[ (\mathrm{m/s})^{-1} \right]$",labelpad=2)

    #plt.subplots_adjust(right=0.97, left=0.15, top=0.95, bottom=0.175)
    plt.subplots_adjust(right=0.97, left=0.13, top=0.92, bottom=0.175)
    plt.savefig(ofoldername + "/vHist_%08d.png" %(outIdx,),dpi=300)
    #plt.show()
    plt.clf()
    
    pCyclic += 1
    outIdx += 1

#Final summary plots
plt.figure(dpi=DPI,figsize=FIGSIZE)
plt.semilogy(log_time, log_Tz, label="$T_z$", color="blue")
plt.semilogy(log_time, log_Tr, label="$T_r$", color="green")
plt.semilogy(log_time, log_Tt, label="$T_\\theta$", color="red")
plt.legend(loc=2)
plt.xlabel("Time [ns]")
plt.ylabel("Temperature [K]")
plt.title("Temperature log for " + species +
          ", region (Z,R) [um]: (%.1f,%.1f)-(%.1f,%.1f)" % (min_z, min_r, max_z, max_r))
plt.savefig(ofoldername + "/summaryplot_temperature.png", dpi=300)
#plt.clf()

plt.figure(dpi=DPI,figsize=FIGSIZE)
plt.semilogy(log_time, log_N)
plt.xlabel("Time [ns]")
plt.ylabel("Number of particles")
plt.title("Number log for " + species +
          ", region (Z,R) [um]: (%.1f,%.1f)-(%.1f,%.1f)" % (min_z, min_r, max_z, max_r))
plt.savefig(ofoldername + "/summaryplot_number_log.png", dpi=300)
#plt.clf()

plt.figure(dpi=DPI,figsize=FIGSIZE)
plt.plot(log_time, log_N)
plt.xlabel("Time [ns]")
plt.ylabel("Number of particles")
plt.title("Number log for " + species +
          ", region (Z,R) [um]: (%.1f,%.1f)-(%.1f,%.1f)" % (min_z, min_r, max_z, max_r))
plt.savefig(ofoldername + "/summaryplot_number_linar.png", dpi=300)
#plt.clf()

plt.figure(dpi=DPI,figsize=FIGSIZE)
plt.semilogy(log_time, log_dens)
plt.xlabel("Time [ns]")
plt.ylabel("Particle number density [cm$^{-3}$]")
plt.title("Density log for " + species +
          ", region (Z,R) [um]: (%.1f,%.1f)-(%.1f,%.1f)" % (min_z, min_r, max_z, max_r))
plt.savefig(ofoldername + "/summaryplot_density_log.png", dpi=300)
#plt.clf()

plt.figure(dpi=DPI,figsize=FIGSIZE)
plt.plot(log_time, log_dens)
plt.xlabel("Time [ns]")
plt.ylabel("Particle number density [cm$^{-3}$]")
plt.title("Density log for " + species +
          ", region (Z,R) [um]: (%.1f,%.1f)-(%.1f,%.1f)" % (min_z, min_r, max_z, max_r))
plt.savefig(ofoldername + "/summaryplot_density_linear.png", dpi=300)
#plt.clf()

#Movie assembly (may crash!)
exit(0) #skip
movieFileName = ofoldername + "/" + os.path.basename(ofoldername) + ".mp4"
if mintime>0:
    out_dt = (timestamp[2]-timestamp[1])*skipFrame
else:
    out_dt = (timestamp[1]-timestamp[0])*skipFrame

fps = speed/out_dt
print "Assembling movie '" + movieFileName + "' at", fps, "fps:"

os.system("rm " + movieFileName)
os.system("ffmpeg -sameq -r "+str(fps)+" -i "+ ofoldername +"/vHist_%08d.png " + movieFileName)
#os.system("rm _tmp*.png")
