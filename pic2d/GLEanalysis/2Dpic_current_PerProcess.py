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
# 2Dpic_current_PerProcess.py:
# Plot current per injection process as a function of time
# Data provided by model 'ArcOriginalNewHS',
# which writes file 'arcbound_original.dat' in the expected format.
#

import sys, os, shutil
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm 
#from matplotlib.ticker import LogFormatter

#from matplotlib import gridspec

from matplotlib import rcParams, rc
rcParams.update({'text.usetex': True})

PLOTOPTS = 1 

if PLOTOPTS == 1: #spread_24x6um_RC_1740V_1pF_1000ohm_bf2_halfgrid
    DPI = 500
    rcParams.update({'savefig.dpi':DPI})
    rc('font',**{'family':'serif','serif':['Times'],'size':8})
    textwidth = 449.40363*1.0/72.27 #inches, for 2Dpic paper
    FIGSIZE = (0.5*0.95*textwidth, 0.5*0.98*textwidth/1.61803398875)
elif PLOTOPTS == 2: #spread_24x6um_RC_1740V_1pF_1000ohm_bf2_halfgrid_Y2
    DPI = 100#500
    rcParams.update({'savefig.dpi':DPI})
    rc('font',**{'family':'serif','serif':['Times'],'size':8})
    textwidth = 449.40363*1.0/72.27 #inches, for 2Dpic paper
    FIGSIZE = None#(0.5*0.95*textwidth, 0.5*0.98*textwidth/1.61803398875)
else:
    DPI = 300
    rcParams.update({'savefig.dpi':DPI})
    FIGSIZE = None

#Get the scaling setup
MODLOAD_parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,MODLOAD_parentdir) 
from calcScaling import InputFile

inputfile = InputFile("../input.txt")
inputfile.calcBasicParams()

if inputfile.ArcBoundaryName != "ArcOriginalNewHS":
    print "ArcBoundaryName must be 'ArcOriginalNewHS' for this analysis to work."


print
print "Usage: 2Dpic_current_PerProcess.py smoothsteps"
smoothsteps = None
if len(sys.argv) == 2:
    smoothsteps = int(sys.argv[1])
    print "- smoothsteps =", smoothsteps
else:
    print "Please specify smoothsteps."
    exit(1)
print

def patchArray(ts, columns,patchmode='zero',dt=1):
    #Used to "patch" up an array with missing data, i.e. if the disk ran full during simulation
    assert dt == 1
    
    fragmentIdx = [0]
    prevTS = ts[0]
    for i in xrange(1,len(ts)):
        if ts[i] != prevTS+dt:
            print "\t Found jump != dt, prevTS="+str(prevTS)+", ts[i]="+str(ts[i])+", i=",i
            fragmentIdx.append(i)
        prevTS = ts[i]
        #fragmentIdx.append(len(i)-1)
    newTS = []
    newColumns = []
    for col in columns:
        newColumns.append([])
    newColumns = tuple(newColumns)
    
    #Patch together fragments and fillers
    print "\t fragmentIdx = " + str(fragmentIdx) + ", len(ts)=" + str(len(ts))
    for i in xrange(1,len(fragmentIdx)):
        print "\t appending, i=" + str(i),
        newTS += list(ts[fragmentIdx[i-1]:fragmentIdx[i]])
        addRange = range(ts[fragmentIdx[i]-1],ts[fragmentIdx[i]],dt)
        newTS += addRange
        print "len(addRange) = " + str(len(addRange))

        for (colOld,colNew) in zip(columns,newColumns):
            colNew += list(colOld[fragmentIdx[i-1]:fragmentIdx[i]])
            if patchmode=='zero':
                colNew += [0]*(len(addRange))
    #tail of arrays
    newTS += list(ts[fragmentIdx[-1]:])
    for (colOld,colNew) in zip(columns,newColumns):
        colNew += list(colOld[fragmentIdx[-1]:])

    return (np.asarray(newTS), map(np.asarray, newColumns) )
    pass

#Read arcbound_original.dat
print "Reading arcbound_original.dat"
(ts, emit_tip, emit_flat, emit_SEY, emit_evap, emit_sput_cat, emit_sput_ano, emit_htspk) = np.loadtxt("../arcbounds_original.dat", dtype='int', converters={2:lambda x: None}, usecols=(0,3,5,6,7,8,9,10), unpack=True)
if PLOTOPTS == 2:
    print "patching data..."
    (ts, (emit_tip, emit_flat, emit_SEY, emit_evap, emit_sput_cat, emit_sput_ano, emit_htspk)) = patchArray(ts, (emit_tip, emit_flat, emit_SEY, emit_evap, emit_sput_cat, emit_sput_ano, emit_htspk))
print "done."
print


#Read circuit.dat
print "Reading circuit.dat"
(circ_ts, circ_deltaQ, circ_U) = np.loadtxt("../circuit.dat", usecols=(0,1,2), unpack=True)
if PLOTOPTS == 2:
    print "patching data..."
    (circ_ts, (circ_deltaQ, circ_U)) = patchArray(circ_ts, (circ_deltaQ, circ_U))
print "done."
print

#Read arcbound.dat
print "Reading arcbound.dat"
(ab_ts, ab_removed_e, ab_removed_i,ab_removed_n_cat,ab_removed_n_ano,ab_removed_n_r) =  np.loadtxt("../arcbounds.dat", usecols=(0,9,15,18,19,20), unpack=True)
if PLOTOPTS == 2:
    print "patching data..."
    (ab_ts, (ab_removed_e, ab_removed_i,ab_removed_n_cat,ab_removed_n_ano,ab_removed_n_r)) = patchArray(ab_ts, (ab_removed_e, ab_removed_i,ab_removed_n_cat,ab_removed_n_ano,ab_removed_n_r))
print "done."
print

#Read mainstats.dat
print "Reading mainStats.dat"
(ms_ts, ms_numNeutrals) = np.loadtxt("../mainStats.dat", usecols=(0,5), unpack=True)
if PLOTOPTS == 2:
    print "patching data..."
    (ms_ts, (ms_numNeutrals,)) = patchArray(ms_ts, (ms_numNeutrals,))
print "done."
print

#Ionization rate
print "Calculating ionization rate..."
assert inputfile.n2inj_step == inputfile.dt_ion
print "dt_ion =", inputfile.dt_ion
delta_numNeutrals = np.diff([0]+list(ms_numNeutrals[::inputfile.dt_ion]))
ionizations = (emit_evap[::inputfile.dt_ion]+emit_sput_cat[::inputfile.dt_ion]+emit_sput_ano[::inputfile.dt_ion]+emit_htspk[::inputfile.dt_ion] - ab_removed_n_cat[::inputfile.dt_ion]-ab_removed_n_ano[::inputfile.dt_ion]-ab_removed_n_r[::inputfile.dt_ion] - delta_numNeutrals)/float(inputfile.dt_ion)
print "done."
print

#Important utility functions
def smooth(data,steps=smoothsteps, endFix=True):
    "Boxcar window smoother"
    #print "Convolving..."
    if smooth.kernel==None or len(smooth.kernel) != steps:
        print "Creating kernel with steps =", \
            steps, "=", steps*inputfile.dT*1e9, "ns at dt,", steps*inputfile.dT*1e9*inputfile.dt_ion, "ns at dt_ion"
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

def backfill(data,steps):
    assert type(data) == np.ndarray #this type passed by reference
    fSteps = float(steps)
    #Normalize first point
    data[0] /= fSteps
    #Backfill
    i = steps
    while i < len(data):
        DATA = data[i]/fSteps
        j = i-steps+1
        for k in xrange(j,i+1):
            data[k] = DATA
        i += steps

print "backfilling..."
dt_ion = inputfile.dt_ion
assert ts[1]-ts[0] == 1 #dt_out
assert inputfile.e2inj_step == 1
assert inputfile.n2inj_step == dt_ion
backfill(emit_evap,     dt_ion)
backfill(emit_sput_cat, dt_ion)
backfill(emit_sput_ano, dt_ion)
backfill(emit_htspk,    dt_ion)
backfill(ab_removed_i,  dt_ion)
backfill(ab_removed_n_cat,  dt_ion)
backfill(ab_removed_n_ano,  dt_ion)
backfill(ab_removed_n_r,  dt_ion)
print "Done."

convFactor1 = inputfile.N_sp/inputfile.dT #superpartices/timestep -> particles/sec
convFactor2 = 1.60217657e-19*inputfile.N_sp/inputfile.dT # superparticles/timestep -> Amps
convFactor3 = (inputfile.dz/inputfile.Omega_pe)**2*inputfile.T_ref #dimless potential -> Volts

print "Smoothing..."
(emit_tip_smooth, emit_flat_smooth, emit_SEY_smooth, \
     emit_evap_smooth, emit_sput_cat_smooth, emit_sput_ano_smooth, emit_htspk_smooth, ab_removed_e_smooth, ab_removed_i_smooth, ab_removed_n_cat_smooth,ab_removed_n_ano_smooth,ab_removed_n_r_smooth) = \
     map(smooth, (emit_tip, emit_flat, emit_SEY, emit_evap, emit_sput_cat, emit_sput_ano, emit_htspk, ab_removed_e, ab_removed_i, ab_removed_n_cat,ab_removed_n_ano,ab_removed_n_r) )
circ_deltaQ_smooth = smooth(circ_deltaQ)
ionizations_smooth = smooth(ionizations,int(1.0*smoothsteps/inputfile.dt_ion))
ionizations_time = ms_ts[::inputfile.dt_ion]*inputfile.dT*1e9
print "Done."


t = ts*inputfile.dT*1e9 #ns
circ_t = circ_ts*inputfile.dT*1e9 #ns
ab_t = ab_ts*inputfile.dT*1e9 #ns

# plt.figure(1) #Particle currents (particles/sec)
# plt.plot(t, emit_tip*convFactor1,       label="Tip")
# plt.plot(t, emit_flat*convFactor1,      label="Flat")
# plt.plot(t, emit_SEY*convFactor1,       label="SEY")
# plt.plot(t, emit_evap*convFactor1,      label="Evap")
# plt.plot(t, emit_sput_cat*convFactor1,  label="Sput (cat)")
# plt.plot(t, emit_sput_ano*convFactor1,  label="Sput (ano)")
# plt.plot(t, emit_htspk*convFactor1,     label="Heatspike")
# plt.plot(circ_t, circ_deltaQ*convFactor1, "r--", label="Circuit current")
# plt.title("Cathode")
# plt.xlabel("t [ns]")
# plt.ylabel("particles/sec.")
# plt.legend(loc=0)

# plt.figure(2,dpi=DPI,figsize=FIGSIZE) #Particle currents (particles/sec) (smoothed)
# plt.plot(t, emit_tip_smooth*convFactor1,       label="Tip")
# plt.plot(t, emit_flat_smooth*convFactor1,      label="Flat")
# plt.plot(t, emit_SEY_smooth*convFactor1,       label="SEY")
# plt.plot(t, emit_evap_smooth*convFactor1,      label="Evap")
# plt.plot(t, emit_sput_cat_smooth*convFactor1,  label="Sput (cat)")
# plt.plot(t, emit_sput_ano_smooth*convFactor1,  label="Sput (ano)")
# plt.plot(t, emit_htspk_smooth*convFactor1,     label="Heatspike")
# plt.plot(circ_t, circ_deltaQ_smooth*convFactor1, "r--", label="Circuit current")
# plt.title("Cathode, smoothing = %d steps = %.3f ps " %
#           (smoothsteps, smoothsteps*inputfile.dT*1e12) )
# plt.xlabel("t [ns]")
# plt.ylabel("particles/sec.")
# plt.legend(loc=0)

plt.figure(3,dpi=DPI,figsize=FIGSIZE) #Neutral particle currents (particles/sec) (smoothed)
line1, = plt.plot(t, emit_evap_smooth*convFactor1,                   label="Evaporation")
line2, = plt.plot(t, emit_sput_cat_smooth*convFactor1,               label="Sputtering from cathode")
line3, = plt.plot(t, emit_sput_ano_smooth*convFactor1,               label="Sputtering from anode")
line4, = plt.plot(t, emit_htspk_smooth*convFactor1,                  label="High-flux sputtering")
#line5, = plt.plot(t, -ab_removed_n_cat_smooth*convFactor1,           label="Removed (cathode)")
#line6, = plt.plot(t, -ab_removed_n_ano_smooth*convFactor1,           label="Removed (anode)")
#line7, = plt.plot(t, -ab_removed_n_r_smooth*convFactor1,             label="Removed (radial)")
if PLOTOPTS == 1:
    #Problematic when missing different ammounts of data in different files
    line8, = plt.plot(ab_t, +(ab_removed_n_cat_smooth+ab_removed_n_ano_smooth+ab_removed_n_r_smooth)*convFactor1,
                      label="Removed at walls \n (no arrow)") 
line9, = plt.plot(ionizations_time, +ionizations_smooth*convFactor1, label="Ionizations")
#
#line10, = plt.plot(circ_t, circ_deltaQ_smooth*convFactor1, "r--", label="Circuit current")
# plt.title("Cathode, smoothing = %d steps = %.3f ps " %
#           (smoothsteps, smoothsteps*inputfile.dT*1e12) )
plt.xlabel("Time [ns]")
plt.ylabel("Particles injected or removed /sec.")
plt.legend(loc=2,frameon=False,fontsize=8)

if PLOTOPTS == 1:
    plt.subplots_adjust(right=0.99, left=0.13, top=0.92, bottom=0.16)

    plt.xlim(0.8,2.1)
    plt.annotate(
        r"Evap.", xy = (1.74,6.3e17), xycoords='data',
        color=line1.get_color(), xytext=(-5, +20), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"),
        horizontalalignment='center', verticalalignment='center',fontsize=8
        )
    plt.annotate(
        r"Hi-f.sput.", xy = (1.74,5.7e17), xycoords='data',
        color=line4.get_color(), xytext=(-5, -20), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"),
        horizontalalignment='center', verticalalignment='center',fontsize=8
        )
    plt.annotate(
        r"Ionizations", xy = (1.79,1.05e18), xycoords='data',
        color=line9.get_color(), xytext=(20, 20), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"),
        horizontalalignment='center', verticalalignment='center',fontsize=8
        )
    plt.annotate(
        r"Sput./ano.", xy = (2.05,8.0e16), xycoords='data',
        color=line3.get_color(), xytext=(-11, 0), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"),
        horizontalalignment='right', verticalalignment='bottom',fontsize=8
        )
    plt.annotate(
        r"Sput./cat.", xy = (1.66,1.4e17), xycoords='data',
        color=line2.get_color(), xytext=(-11, 0), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"),
        horizontalalignment='right', verticalalignment='center',fontsize=8
        )

plt.savefig("pngs/current_PerProcess_neutralFluxes_smoothed.png",dpi=DPI)

# plt.figure(4) #Charged particle currents (amps)
# plt.plot(t, emit_tip*convFactor2,       label="Tip")
# plt.plot(t, emit_flat*convFactor2,      label="Flat")
# plt.plot(t, emit_SEY*convFactor2,       label="SEY")
# plt.plot(circ_t, circ_deltaQ*convFactor2, "--", label="Circuit current")
# plt.title("Cathode")
# plt.xlabel("t [ns]")
# plt.ylabel("Current [A]")
# plt.legend(loc=0)

plt.figure(5, dpi=DPI,figsize=FIGSIZE) #Particle currents (amps) (smoothed)
line1, = plt.plot(t, emit_tip_smooth*convFactor2,       label="Tip")
line2, = plt.plot(t, emit_flat_smooth*convFactor2,      label="Flat")
line3, = plt.plot(t, emit_SEY_smooth*convFactor2,       label="SEY")
line4, = plt.plot(ab_t, ab_removed_i_smooth*convFactor2,   label=r"Removed Cu\textsuperscript{+}")
line5, = plt.plot(ab_t, -1*ab_removed_e_smooth*convFactor2,label=r"Removed e\textsuperscript{-}")
line6, = plt.plot(circ_t, circ_deltaQ_smooth*convFactor2, ":", label=r"$I_\mathrm{circ}$")
# plt.title("Cathode, smoothing = %d steps = %.3f ps " %
#           (smoothsteps, smoothsteps*inputfile.dT*1e12) )
plt.xlabel("Time [ns]")
plt.ylabel("Current [A]")
plt.legend(loc=2,frameon=False,ncol=1,fontsize=8)

if PLOTOPTS == 1:
    #plt.subplots_adjust(right=0.99, left=0.13, top=0.97, bottom=0.16)
    plt.subplots_adjust(right=0.88, left=0.13, top=0.92, bottom=0.16) #make room for extra y-axis
    
    plt.xlim(0.8,2.1)
    plt.annotate(
        r"Flat", xy = (1.4,1.9), xycoords='data',
        color=line2.get_color(), xytext=(0, -12), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"),
        horizontalalignment='center', verticalalignment='center',fontsize=8
        )
    plt.annotate(
        r"$I_\mathrm{circ}$", xy = (1.47,2.5), xycoords='data',
        color=line6.get_color(), xytext=(0, +20), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"),
        horizontalalignment='center', verticalalignment='center',fontsize=8
        )
    plt.annotate(
        r"Tip", xy = (1.93,0.33), xycoords='data',
        color=line1.get_color(), xytext=(-10, +20), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"),
        horizontalalignment='right', verticalalignment='center',fontsize=8
        )
    plt.annotate(
        r"Removed e\textsuperscript{-}", xy = (1.62,-0.5), xycoords='data',
        color=line5.get_color(), xytext=(-10, -10), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"),
        horizontalalignment='right', verticalalignment='center',fontsize=8
        )
    plt.annotate(
        r"Removed Cu\textsuperscript{+}", xy = (1.65,0.0), xycoords='data',
        color=line4.get_color(), xytext=(0, 17), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"),
        horizontalalignment='center', verticalalignment='center',fontsize=8
        )
    plt.annotate(
        r"SEY", xy = (1.75,0.0), xycoords='data',
        color=line3.get_color(), xytext=(10, -18), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"),
        horizontalalignment='left', verticalalignment='center',fontsize=8
        )

    #right-hand y-axis with alternative units
    print "Adding right hand y-axis"
    (xmin,xmax,ymin,ymax) = plt.axis()
    ax1=plt.gca()
    ax1.yaxis.tick_left()
    ax2 = ax1.twinx()
    plt.ylabel("Particles injected or removed/sec.")
    ax2.yaxis.tick_right()
    plt.axis((xmin,xmax,ymin/1.60217657e-19,ymax/1.60217657e-19))

plt.savefig("pngs/current_PerProcess_cathodeAmps_smoothed.png",dpi=DPI)

# plt.figure(6) #Current and voltage
# plt.gca().plot(circ_t, circ_deltaQ*convFactor2, "b-", label="Circuit current")
# plt.gca().set_ylabel("Current [A]", color='b')
# plt.gca().set_xlabel("t [ns]")
# for tl in plt.gca().get_yticklabels():
#     tl.set_color('b')
# ax2 = plt.gca().twinx()
# ax2.plot(circ_t, circ_U*convFactor3, "r-", label="Circuit current")
# ax2.set_ylabel("Voltage [V]", color='r')
# for tl in ax2.get_yticklabels():
#     tl.set_color('r')


#plt.show()
