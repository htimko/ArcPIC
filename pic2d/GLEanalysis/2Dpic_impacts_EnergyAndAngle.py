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
# 2Dpic_impacts_EnergyAndAngle.py:
# Plots the energy and incidence angle distribution of impacting particles
#


import sys, os, shutil
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm 
from matplotlib.ticker import MultipleLocator

#from matplotlib import gridspec

from matplotlib import rcParams,rc
rcParams.update({'text.usetex': True})
DPI = 500
rcParams.update({'savefig.dpi':DPI})
rc('font',**{'family':'serif','serif':['Times'],'size':8})
textwidth = 449.40363*1.0/72.27 #inches, for 2Dpic paper
FIGSIZE = (0.5*0.95*textwidth, 0.5*0.98*textwidth/1.61803398875)



#Get input
if len(sys.argv) != 7:
    print "Usage: ./2Dpic_impacts_EnergyAndAngle.py <e|i|n> <c|a|r> <mintime> <maxtime> nbinsX nbinsY"
    exit(1)

species = sys.argv[1]
if not (species == "e" or species == "i" or species == "n"):
    print "Species must be 'e', 'i', or 'n'"
    exit(1)

edg = sys.argv[2]
if not (edg == "c" or edg == "a" or edg == "r"):
    print "Edge must be 'c', 'a', or 'r'"
    exit(1)
if (edg == 'r'):
    print "Edge 'r' not supported at this time"
    exit(1)

mintime = int(sys.argv[3])
maxtime = int(sys.argv[4])
nbinsX  = int(sys.argv[5])
nbinsY  = int(sys.argv[6])

# minR = None
# maxR = None

print "Got options:"
print " - mintime = ", mintime
print " - maxtime = ", maxtime
print
print " - species = ", species
print " - edge    = ", edg
print
print " - nbinsX  = ", nbinsX
print " - nbinsY  = ", nbinsY
print

#Get the scaling setup
MODLOAD_parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,MODLOAD_parentdir) 
from calcScaling import InputFile
inputfile = InputFile("../input.txt")
inputfile.calcBasicParams()
print
print "Time window =", mintime*inputfile.dT*1e9, "--", maxtime*inputfile.dT*1e9, "ns"
print

# if minR == maxR and minR == None:
#     minR = 0.0
#     maxR = inputfile.nr

step = []
z    = []
r    = []
vz   = []
vr   = []
vt   = []
edge = [] #c = cathode, a = anode, r = radial

maxZ = float(inputfile.nz)

particleFileName = "removed_" + species + ".dat"
print "Reading file '%s'..." % (particleFileName,)
totalcount = 0
particleFile = open("../out/" + particleFileName, 'r')
particleFile.readline() #Skip comments line
for line in particleFile:
    l = line.split()
    if len(l) == 0:
        print "DROP empty '%s'" %(line,)
        continue
    STEP = int(l[0])
    Z    = float(l[1])
    R    = float(l[2])
    VZ   = float(l[3])
    VR   = float(l[4])
    VT   = float(l[5])

    if STEP < mintime:
        #print "DROP step", line
        continue
    if STEP > maxtime:
        #print "STOP step", line
        break

    totalcount += 1
    if totalcount % 1000000 == 0:
        print "Read", totalcount/1000000, "million particles, step=", STEP, "=", STEP*inputfile.dT*1e9, "[ns]"

    # if R < minR or R > maxR:
    #     print "DROP R", line
    #     continue

    if Z <= 0.0:
        EDGE = "c"
    elif Z >= maxZ:
        EDGE = "a"
    else:
        EDGE = "r"
    if EDGE != edg:
        #print "DROP EDG", line,
        continue

#    print "KEEP", line,

    step.append(STEP)
    z.append(Z)
    r.append(R)
    vz.append(VZ)
    vr.append(VR)
    vt.append(VT)
    edge.append(EDGE);
    
particleFile.close()
print "done."

if len(step) == 0:
    print "No data collected, check window."
    exit(1)
print "Got", len(step), "particles, rejected", totalcount-len(step)
print

print "Postprocessing..."

#Convert
posNorm = inputfile.dZ*1e4 # convert to um
velNorm = inputfile.dZ*1e-2 / inputfile.dT # convert to m/s
if species != 'e':
    velNorm /= float(inputfile.dt_ion)

mass = 510.998928e3 #electron mass, eV/c^2 
if species == 'n':
    mass =  63.546/1.0073 * inputfile.mi_over_me*mass
elif species == 'i':
    mass = (63.546/1.0073 * inputfile.mi_over_me-1)*mass

z = np.asarray(z)*posNorm
r = np.asarray(r)*posNorm

vz = np.asarray(vz)*velNorm
vr = np.asarray(vr)*velNorm
vt = np.asarray(vt)*velNorm

step = np.asarray(step)
edge = np.asarray(edge)

time = step*inputfile.dT*1e9 #[ns]

if maxtime > step[-1]:
    maxtime = step[-1]

angle = np.arctan(vr/vz)*180.0/np.pi
if edg == 'c':
    angle *= -1.0

energy = 0.5*mass/299792458.0**2 * (vz**2 + vr**2 + vt**2) #eV

if species != 'e':
    print "Calculating Y-T yield...",
    YT_yield_eps = (4.45394e-06)*energy;
    YT_yield_sn  = 8205*3.441*np.sqrt(YT_yield_eps)*np.log(YT_yield_eps + 2.718) / (1 + 6.355*np.sqrt(YT_yield_eps) + YT_yield_eps*(6.882*np.sqrt(YT_yield_eps) - 1.708));
    YT_yield = 0.042*(0.2525/3.49) * (YT_yield_sn/(1 + 0.010124*0.1573*np.power(YT_yield_eps,0.3))) * np.power(1 - np.sqrt(23.383/energy),2.5);
    for i in xrange(len(YT_yield)):
        if energy[i] < 23.383:
            YT_yield[i] = 0.0
    del YT_yield_eps
    del YT_yield_sn
    print "done."
    print "Average yield = ", np.average(YT_yield)

print "Plotting:"

# print "Figure 1"
# plt.figure(1,dpi=DPI,figsize=FIGSIZE)
# plt.hist(angle, nbinsX*nbinsY, lw=0)
# plt.xlabel("Incidence angle [degrees]")
# plt.ylabel("Counts")
# #plt.title("Species: %c, surface: %c, time window: %f-%f ns" %(species, edg, mintime*inputfile.dT*1e9, maxtime*inputfile.dT*1e9) )

print "Figure 2"
plt.figure(2,dpi=DPI,figsize=FIGSIZE)
histweights=np.ones(len(energy))/((np.max(energy)-np.min(energy))/(nbinsX*nbinsY))
plt.hist(energy, nbinsX*nbinsY, lw=0, range=(np.min(energy),np.max(energy)), weights=histweights )
plt.xlabel("Incidence energy [eV]")
plt.ylabel("Number of particles / eV")
#plt.title("Species: %c, surface: %c, time window: %f-%f ns" %(species, edg, mintime*inputfile.dT*1e9, maxtime*inputfile.dT*1e9) )
plt.xlim(0,205)
plt.gca().xaxis.set_major_locator(MultipleLocator(25))
plt.gca().xaxis.set_minor_locator(MultipleLocator(5))
plt.subplots_adjust(right=0.99, left=0.2, top=0.95, bottom=0.17)

# print "Figure 3"
# plt.figure(3,dpi=DPI,figsize=FIGSIZE)
# plt.hist2d(energy,r, bins=(nbinsX, nbinsY), 
#            range=((0,max(energy)),(0,inputfile.R*1e4)), norm=LogNorm())
# plt.colorbar().set_label("particles / bin")
# plt.xlabel("Incidence energy [eV]")
# plt.ylabel("Radius [um]")
# #plt.title("Species: %c, surface: %c, time window: %f-%f ns" %(species, edg, mintime*inputfile.dT*1e9, maxtime*inputfile.dT*1e9) )

# print "Figure 4"
# plt.figure(4,dpi=DPI,figsize=FIGSIZE)
# plt.hist2d(energy, time, bins=(nbinsX,nbinsY),norm=LogNorm())
# plt.colorbar().set_label("particles/bin")
# plt.xlabel("Incidence energy [eV]")
# plt.ylabel("Time [ns]")
# #plt.title("Species: %c, surface: %c, time window: %f-%f ns" %(species, edg, mintime*inputfile.dT*1e9, maxtime*inputfile.dT*1e9) )

# print "Figure 5"
# plt.figure(5,dpi=DPI,figsize=FIGSIZE) #same as 4 but flip x/y
# plt.hist2d(time, energy, bins=(nbinsX,nbinsY),norm=LogNorm())
# plt.colorbar().set_label("particles / bin")
# plt.xlabel("Time [ns]")
# plt.ylabel("Incidence energy [eV]")
# #plt.title("Species: %c, surface: %c, time window: %f-%f ns" %(species, edg, mintime*inputfile.dT*1e9, maxtime*inputfile.dT*1e9) )

print "Figure 6"
plt.figure(6,dpi=DPI,figsize=FIGSIZE)
(trCounts, trXedges, trYedges, trImage) = \
    plt.hist2d(time, r, bins=(nbinsX,nbinsY),
               range=((0.0,time[-1]),(0.0,inputfile.R*1e4)),norm=LogNorm())
plt.colorbar().set_label("particles / bin")
plt.xlabel("Time [ns]")
plt.ylabel("Radius [um]")
#plt.title("Species: %c, surface: %c, time window: %f-%f ns" %(species, edg, mintime*inputfile.dT*1e9, maxtime*inputfile.dT*1e9) )

print "Figure 7"
plt.figure(7,dpi=DPI,figsize=FIGSIZE)
dt = (trXedges[-1]-trXedges[0])/(len(trXedges)-1.0)*1e-9 #s
trCounts = trCounts.copy()
#print trCounts
for i in xrange(len(trYedges)-1):
    r1 = trYedges[i]*1e-4 #cm
    r2 = trYedges[i+1]*1e-4   #cm
    area = np.pi*(r2**2-r1**2) #cm^2
    #print r1, r2, area
    trCounts[:,i] /= area*dt #particles/s/cm^2
#print trCounts
plt.pcolormesh(trXedges,trYedges,trCounts.transpose(), norm=LogNorm())
plt.ylim(0.0,inputfile.R*1e4)
plt.colorbar().set_label("Flux density [Particles / s / cm$^2$]")
plt.xlabel("Time [ns]")
plt.ylabel("Radius [\\textmu m]")
#plt.title("Species: %c, surface: %c, time window: %f-%f ns" %(species, edg, mintime*inputfile.dT*1e9, maxtime*inputfile.dT*1e9) )
plt.xlim(0,2.1)
plt.subplots_adjust(right=0.95, left=0.12, top=0.95, bottom=0.17)

# print "Figure 8"
# plt.figure(8,dpi=DPI,figsize=FIGSIZE)
# currentAdd = np.empty(len(trYedges)-1)
# for i in xrange(len(trYedges)-1):
#     currentAdd[i] = np.sum(trCounts[:,i])
#     currentAdd[i] /= len(trXedges)-1.0
# plt.bar(trYedges[:-1], currentAdd, width=(trYedges[-1]-trYedges[0])/(len(trYedges)-1.0), lw=0.0, log=True)
# plt.xlabel("Radius[um]")
# plt.ylabel("Flux [Particles / s / cm$^2$]");
# #plt.title("Species: %c, surface: %c, time window: %f-%f ns" %(species, edg, mintime*inputfile.dT*1e9, maxtime*inputfile.dT*1e9) )

# print "Figure 9"
# plt.figure(9,dpi=DPI,figsize=FIGSIZE)
# plt.hist(YT_yield, nbinsX*nbinsY, lw=0)
# plt.axvline(1.0, ls='--', color='k')
# plt.xlabel("Sputtering yield");
# plt.ylabel("Counts")
# #plt.title("Species: %c, surface: %c, time window: %f-%f ns" %(species, edg, mintime*inputfile.dT*1e9, maxtime*inputfile.dT*1e9) )

# print "Figure 10"
# plt.figure(10,dpi=DPI,figsize=FIGSIZE)
# plt.hist2d(time, YT_yield, bins=(nbinsX,nbinsY),norm=LogNorm())
# plt.colorbar().set_label("particles / bin")
# plt.axhline(1.0, ls='--', color='k')
# plt.xlabel("Time [ns]")
# plt.ylabel("Sputtering yield")
# #plt.title("Species: %c, surface: %c, time window: %f-%f ns" %(species, edg, mintime*inputfile.dT*1e9, maxtime*inputfile.dT*1e9) )


plt.show()
