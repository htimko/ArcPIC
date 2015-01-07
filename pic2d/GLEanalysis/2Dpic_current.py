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
# 2Dpic_current.py:
# Plot current as a function of time
#

import sys, os, shutil
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm 
#from matplotlib.ticker import LogFormatter

#from matplotlib import gridspec

from matplotlib import rcParams
rcParams.update({'text.usetex': True})

#Get the scaling setup
MODLOAD_parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,MODLOAD_parentdir) 
from calcScaling import InputFile

inputfile = InputFile("../input.txt")
inputfile.calcBasicParams()

print
print "Usage: 2Dpic_current.py smoothsteps (maxsteps)"
smoothsteps = None
maxsteps = None
if len(sys.argv) == 2 or len(sys.argv) == 3:
    smoothsteps = int(sys.argv[1])
    print "- smoothsteps =", smoothsteps
else:
    print "Please specify smoothsteps."
    exit(1)
if len(sys.argv) == 3:
    maxsteps = int(sys.argv[2])
    print "- maxsteps    =", maxsteps, "=", maxsteps*inputfile.dT*1e9, "ns"
print

#Read arcbounds.dat
ts = []

inj_e_c = []
inj_e_a = []

inj_n_c = []
inj_n_a = []

removed_e_c = []
removed_e_a = []
removed_e_r = []

removed_i_c = []
removed_i_a = []
removed_i_r = []

removed_n_c = []
removed_n_a = []
removed_n_r = []

arcbounds = open("../arcbounds.dat")
arcbounds.readline() #skip header lines
arcbounds.readline()
arcbounds.readline()
arcbounds.readline()

print "parsing arcbounds.dat"
for line in arcbounds:
    l = line.split()
    if len(l) == 0:
        print "Empty line? line='"+line+"'"
        continue

    ts.append(int(l[0]))

    inj_e_c.append(int(l[1]))
    inj_e_a.append(int(l[2]))
    
    inj_n_c.append(int(l[7]))
    inj_n_a.append(int(l[8]))
    
    removed_e_c.append(int(l[9]))
    removed_e_a.append(int(l[10]))
    removed_e_r.append(int(l[11]))
    
    removed_i_c.append(int(l[15]))
    removed_i_a.append(int(l[16]))
    removed_i_r.append(int(l[17]))
    
    removed_n_c.append(int(l[18]))
    removed_n_a.append(int(l[19]))
    removed_n_r.append(int(l[20]))

    if maxsteps != None and ts[-1] >= maxsteps:
        print "maxsteps reached"
        break
arcbounds.close()
print "done."

ts = np.asarray(ts)

inj_e_c = np.asarray(inj_e_c, 'float')
inj_e_a = np.asarray(inj_e_a, 'float')

inj_n_c = np.asarray(inj_n_c, 'float')
inj_n_a = np.asarray(inj_n_a, 'float')

removed_e_c = np.asarray(removed_e_c, 'float')
removed_e_a = np.asarray(removed_e_a, 'float')
removed_e_r = np.asarray(removed_e_r, 'float')

removed_i_c = np.asarray(removed_i_c, 'float')
removed_i_a = np.asarray(removed_i_a, 'float')
removed_i_r = np.asarray(removed_i_r, 'float')

removed_n_c = np.asarray(removed_n_c, 'float')
removed_n_a = np.asarray(removed_n_a, 'float')
removed_n_r = np.asarray(removed_n_r, 'float')

#Smooth the neutral injection
print "neutral injection/removal backfilling"
dt_ion = inputfile.dt_ion
dt_out = ts[1]-ts[0]
assert dt_out == 1
assert inputfile.e2inj_step == 1
assert inputfile.n2inj_step == dt_ion
if dt_ion == 1:
    print "dt_ion = 1, skipping backfill."
else:
    #Normalize first point
    inj_n_c[0] /= float(dt_ion)
    inj_n_a[0] /= float(dt_ion)

    removed_i_c[0] /= float(dt_ion)
    removed_i_a[0] /= float(dt_ion)
    removed_i_r[0] /= float(dt_ion)

    removed_n_c[0] /= float(dt_ion)
    removed_n_a[0] /= float(dt_ion)
    removed_n_r[0] /= float(dt_ion)
    
    #Backfill
    i = dt_ion
    while i < len(ts):
        INJ_N_C = inj_n_c[i]/float(dt_ion)
        INJ_N_A = inj_n_a[i]/float(dt_ion)

        REMOVED_I_C = removed_i_c[i]/float(dt_ion)
        REMOVED_I_A = removed_i_a[i]/float(dt_ion)
        REMOVED_I_R = removed_i_r[i]/float(dt_ion)

        REMOVED_N_C = removed_n_c[i]/float(dt_ion)
        REMOVED_N_A = removed_n_a[i]/float(dt_ion)
        REMOVED_N_R = removed_n_r[i]/float(dt_ion)
        
        j = i-dt_ion+1
        for k in xrange(j,i+1):
            inj_n_c[k] = INJ_N_C
            inj_n_a[k] = INJ_N_A
            
            removed_i_c[k] = REMOVED_I_C
            removed_i_a[k] = REMOVED_I_A
            removed_i_r[k] = REMOVED_I_R
            
            removed_n_c[k] = REMOVED_N_C
            removed_n_a[k] = REMOVED_N_A
            removed_n_r[k] = REMOVED_N_R
            
        i += dt_ion
    print "backfiller finished"

def smooth(data,steps=smoothsteps, endFix=True):
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

#Testplot for smoother
# X = np.linspace(0,1,40)
# Y = np.ones(len(X))
# plt.plot(X,Y, '-+')
# plt.plot(X,smooth(Y,smoothsteps,False), '-+')
# plt.plot(X,smooth(Y,smoothsteps), '-*')
# plt.show()
# exit(0)

#Read circuit.dat
print "Reading circuit.dat"
(circ_ts, circ_deltaQ) = np.loadtxt("../circuit.dat", usecols=(0,1), unpack=True)
print "done."
print

convFactor1 = inputfile.N_sp/inputfile.dT #superpartices/timestep -> particles/sec
convFactor2 = 1.60217657e-19*inputfile.N_sp/inputfile.dT # superparticles/timestep -> Amps

I_cathode = convFactor2 * (inj_e_c + removed_i_c - removed_e_c)
I_anode   = convFactor2 * (removed_e_a - inj_e_a - removed_i_a)
I_rad     = convFactor2 * (removed_e_r - removed_i_r)

t = ts*inputfile.dT*1e9 #ns
circ_t = circ_ts*inputfile.dT*1e9 #ns

plt.figure(1) #Wall currents
plt.plot(t,I_cathode, label="Cathode")
plt.plot(t,I_anode,   label="Anode")
plt.plot(t,I_rad,     label="Radial")
plt.plot(circ_t,convFactor2*circ_deltaQ, label="Cathode (S.R.)")
plt.title("Wall currents, no smoothing")
plt.xlabel("Time [ns]")
plt.ylabel("Current [A]")
plt.legend(loc=0)

plt.figure(2) #Wall currents (smoothed)
plt.plot(t,smooth(I_cathode), label="Cathode")
plt.plot(t,smooth(I_anode),   label="Anode")
plt.plot(t,smooth(I_rad),     label="Radial")
plt.title("Wall currents, smoothing = %d steps = %.3f ps " %
          (smoothsteps, smoothsteps*inputfile.dT*1e12) )
plt.plot(circ_t,convFactor2*smooth(circ_deltaQ), label="Cathode (S.R.)")
plt.xlabel("Time [ns]")
plt.ylabel("Current [A]")
plt.legend(loc=0)

#Neutral particle currents (particles / s)
Ip_cathode_n = convFactor1*( inj_n_c-removed_n_c)
Ip_anode_n   = convFactor1*(-inj_n_a+removed_n_a)
Ip_rad_n     = convFactor1*removed_n_r

plt.figure(3) #Neutral currents
plt.plot(t,Ip_cathode_n, label="Cathode")
plt.plot(t,Ip_anode_n, label="Anode")
plt.plot(t,Ip_rad_n, label="Radial")
plt.title("Neutral wall currents, no smoothing\n(positive: cathode to anode \& r)")
plt.xlabel("Time [ns]")
plt.ylabel("particles/sec")
plt.legend(loc=0)

plt.figure(4) #Neutral currents (smoothed)
plt.plot(t,smooth(Ip_cathode_n), label="Cathode")
plt.plot(t,smooth(Ip_anode_n), label="Anode")
plt.plot(t,smooth(Ip_rad_n), label="Radial")
plt.title("Neutral wall currents, smoothing = %d steps = %.3f ps\n" %
          (smoothsteps, smoothsteps*inputfile.dT*1e12) +
          "(positive: cathode to anode \& r)" )
plt.xlabel("Time [ns]")
plt.ylabel("particles/sec")
plt.legend(loc=0)

plt.figure(5) #Cathode current by species and injected/removed
plt.plot(t,inj_e_c*convFactor1, label="Injected electrons")
plt.plot(t,removed_e_c*convFactor1, label="Removed electrons")
plt.plot(t,inj_n_c*convFactor1, label="Injected neutrals")
plt.plot(t,removed_n_c*convFactor1, label="Removed neutrals")
plt.plot(t,removed_i_c*convFactor1, label="Removed ions")
plt.title("Particle wall currents on cathode, no smoothing")
plt.xlabel("Time [ns]")
plt.ylabel("particles/sec")
plt.legend(loc=0)

plt.figure(6) #Cathode current by species and injected/removed, with smoothing
plt.plot(t,smooth(inj_e_c*convFactor1),
         label="Injected electrons")
plt.plot(t,smooth(removed_e_c*convFactor1),
         label="Removed electrons")
plt.plot(t,smooth(inj_n_c*convFactor1),
         label="Injected neutrals")
plt.plot(t,smooth(removed_n_c*convFactor1),
         label="Removed neutrals")
plt.plot(t,smooth(removed_i_c*convFactor1),
         label="Removed ions")
plt.title("Particle wall currents on cathode, smoothing = %d steps = %.3f ps " %
          (smoothsteps, smoothsteps*inputfile.dT*1e12) )
plt.xlabel("Time [ns]")
plt.ylabel("particles/sec")
plt.legend(loc=0)

plt.figure(7) #Anode current by species and injected/removed
plt.plot(t,inj_e_a*convFactor1, label="Injected electrons")
plt.plot(t,removed_e_a*convFactor1, label="Removed electrons")
plt.plot(t,inj_n_a*convFactor1, label="Injected neutrals")
plt.plot(t,removed_n_a*convFactor1, label="Removed neutrals")
plt.plot(t,removed_i_a*convFactor1, label="Removed ions")
plt.title("Particle wall currents on anode, no smoothing")
plt.xlabel("Time [ns]")
plt.ylabel("particles/sec")
plt.legend(loc=0)

plt.figure(8) #Anode current by species and injected/removed, with smoothing
plt.plot(t,smooth(inj_e_a*convFactor1),
         label="Injected electrons")
plt.plot(t,smooth(removed_e_a*convFactor1),
         label="Removed electrons")
plt.plot(t,smooth(inj_n_a*convFactor1),
         label="Injected neutrals")
plt.plot(t,smooth(removed_n_a*convFactor1),
         label="Removed neutrals")
plt.plot(t,smooth(removed_i_a*convFactor1),
         label="Removed ions")
plt.title("Particle wall currents on anode, smoothing = %d steps = %.3f ps " %
          (smoothsteps, smoothsteps*inputfile.dT*1e12) )
plt.xlabel("Time [ns]")
plt.ylabel("particles/sec")
plt.legend(loc=0)

plt.figure(9) #Radial current by species
plt.plot(t,removed_e_r*convFactor1, label="Removed electrons")
plt.plot(t,removed_n_r*convFactor1, label="Removed neutrals")
plt.plot(t,removed_i_r*convFactor1, label="Removed ions")
plt.title("Particle wall currents on radial boundary, no smoothing")
plt.xlabel("Time [ns]")
plt.ylabel("particles/sec")
plt.legend(loc=0)

plt.figure(10) #Radial current by species, with smoothing
plt.plot(t,smooth(removed_e_r*convFactor1),
         label="Removed electrons")
plt.plot(t,smooth(removed_n_r*convFactor1),
         label="Removed neutrals")
plt.plot(t,smooth(removed_i_r*convFactor1),
         label="Removed ions")
plt.title("Particle wall currents on radial boundary, smoothing = %d steps = %.3f ps " %
          (smoothsteps, smoothsteps*inputfile.dT*1e12) )
plt.xlabel("Time [ns]")
plt.ylabel("particles/sec")
plt.legend(loc=0)

plt.show()
