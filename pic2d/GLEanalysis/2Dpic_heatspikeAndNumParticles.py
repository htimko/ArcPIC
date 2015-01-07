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
# 2Dpic_heatspikeAndNumParticles.py:
# Plot number of neutrals together with heatspike data
#

import sys, os, shutil
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm 
#from matplotlib.ticker import LogFormatter

#from matplotlib import gridspec

from matplotlib import rcParams
rcParams.update({'text.usetex': True})

print "Optional argument: Max number of lines to read"
maxLines = None
if len(sys.argv) == 2:
    maxLines = int(sys.argv[1])
    print "Got maxLines =", maxLines
print

#Get the scaling setup
MODLOAD_parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,MODLOAD_parentdir) 
from calcScaling import InputFile
inputfile = InputFile("../input.txt")
inputfile.calcBasicParams()

ms_sn = []
ms_st = []
ms_ne = []
ms_ni = []
ms_nn = []

print "reading mainStats.dat..."
mainstats = open("../mainStats.dat", 'r')
mainstats.readline()
i = 0
for line in mainstats:
    l = line.split()
    ms_sn.append(int(l[0]))
    ms_st.append(float(l[1]))
    ms_ne.append(int(l[3]))
    ms_ni.append(int(l[4]))
    ms_nn.append(int(l[5]))
    
    if maxLines != None and i >= maxLines:
        break
    i += 1
mainstats.close()
ms_sn = np.asarray(ms_sn)
ms_st = np.asarray(ms_st)
ms_ne = np.asarray(ms_ne)
ms_ni = np.asarray(ms_ni)
ms_nn = np.asarray(ms_nn)
print "done."
print

ab_sn = []
ab_beta = []
ab_emitted = []
ab_Etip = []
ab_sigma = []
ab_HSsteps = []

print "reading arcbounds_original.dat.."
arcbounds = open("../arcbounds_original.dat", 'r')
arcbounds.readline()
i = 0
for line in arcbounds:
    l = line.split()
    ab_sn.append(int(l[0]))
    ab_beta.append(float(l[1]))
    ab_emitted.append(int(l[3]))
    ab_Etip.append(float(l[4]))

    if len(l) == 5:
        #old style arcbounds_original.dat
        ab_sigma.append(int(float(l[5])))
    elif len(l) == 13:
        ab_sigma.append(int(float(l[11])))
    else:
        print "unknown arcbounds_original.dat format"
        exit(1);
    if not (ab_sigma[-1] <= 0 or ab_sigma[-1] > 10000):
        ab_HSsteps.append(i)

    if maxLines != None and i >= maxLines:
        break
    i += 1
arcbounds.close()
print "done."
print

ab_sn = np.asarray(ab_sn)
ab_beta = np.asarray(ab_beta)
ab_emitted = np.asarray(ab_emitted)
ab_Etip = np.asarray(ab_Etip)
ab_sigma = np.asarray(ab_sigma)

theLength = None
print "len(ab_sn) =", len(ab_sn), " len(ms_sn) =", len(ms_sn)
if (len(ab_sn) < len(ms_sn)):
    theLength = len(ab_sn)
    print "ab_sn smallest"
elif (len(ab_sn) > len(ms_sn)):
    theLength = len(ms_sn)
    print "ms_sn smallest"
else:
    print "equal length!"
if theLength != None:
    print "cutting to", theLength
    
    ms_sn = ms_sn[:theLength]
    ms_st = ms_st[:theLength]
    ms_ne = ms_ne[:theLength]
    ms_ni = ms_ni[:theLength]
    ms_nn = ms_nn[:theLength]

    ab_sn      = ab_sn[:theLength]
    ab_beta    = ab_beta[:theLength]
    ab_emitted = ab_emitted[:theLength]
    ab_Etip    = ab_Etip[:theLength]
    ab_sigma   = ab_sigma[:theLength]
    
    ab_HSsteps_temp = []
    for i in ab_HSsteps:
        if i < theLength:
            ab_HSsteps_temp.append(i)

if (ms_sn[-1] != ab_sn[-1]):
    print ms_sn[-1], ab_sn[-1]
    print len(ms_sn), len(ab_sn)
    exit(1)

#print "Heatspike indices:"
#print ab_HSsteps

def smooth(data,steps=3000, endFix=True):
    "Boxcar window smoother"
    #print "Convolving..."
    if smooth.kernel==None or len(smooth.kernel) != steps:
        print "Creating kernel with steps =", \
            steps, "=", steps*inputfile.dT*1e9, "ns"
        smooth.kernel = np.ones(steps)/float(steps)
    if len(data) < steps:
        print "Error: len(data) < steps!"
        exit(1)
    ret = np.convolve(data,smooth.kernel,'same')
    if endFix:
        #print ret
        overhang1 = int(steps)/2
        overhang2 = steps-overhang1-1
        #print overhang1, overhang2
        for i in xrange(overhang1):
            ret[i]    *= float(steps) / float( steps - overhang1 + i )
        for i in xrange(overhang2):
            ret[-i-1] *= float(steps) / float( steps - overhang2 + i )
        #print ret
    #print "done."
    return ret
smooth.kernel=None


plt.figure(1)
for i in ab_HSsteps:
    plt.axvline(x=ms_st[i], c='k')
plt.plot(ms_st, ms_ne*inputfile.N_sp, label="Electrons")
plt.plot(ms_st, ms_ni*inputfile.N_sp, label="Ions")
plt.plot(ms_st, ms_nn*inputfile.N_sp, label="Neutrals")
plt.xlabel("t [ns]")
plt.ylabel("Number of particles")
plt.legend(loc=2)

plt.figure(2)
for i in ab_HSsteps:
    plt.axvline(x=ms_st[i], c='k')
plt.plot(ms_st, smooth(ms_ne*inputfile.N_sp), label="Electrons")
plt.plot(ms_st, smooth(ms_ni*inputfile.N_sp), label="Ions")
plt.plot(ms_st, smooth(ms_nn*inputfile.N_sp), label="Neutrals")
plt.xlabel("t [ns]")
plt.ylabel("Number of particles")
plt.legend(loc=2)


plt.figure(3)
for i in ab_HSsteps:
    plt.axvline(x=ms_st[i], c='k')
plt.plot(ms_st, ab_Etip)
plt.xlabel("t [ns]")
plt.ylabel("$\mathrm{E_{tip}}$ [GV/m]")

plt.show()

