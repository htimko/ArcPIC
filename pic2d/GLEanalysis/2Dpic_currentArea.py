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
# 2Dpic_currentArea.py:
# Plots current as a function of area,
# which is estimated from heatspike footing
# or sheeth extent (as estimated by r where field > threshold)
#

import os, sys

import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams.update({'text.usetex': True})

#Get the scaling setup
MODLOAD_parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,MODLOAD_parentdir) 
from calcScaling import InputFile
inputfile = InputFile("../input.txt")
inputfile.calcBasicParams()

convFactor1 = inputfile.N_sp/inputfile.dT #superpartices/timestep -> particles/sec
convFactor2 = 1.60217657e-19*inputfile.N_sp/inputfile.dT # superparticles/timestep -> Amps
convFactor3 = (inputfile.dz/inputfile.Omega_pe)**2*inputfile.T_ref #dimless potential -> Volts

# Load circuit data
print "Loading circuit.dat..."
(circ_ts, circ_deltaQ, circ_U) = np.loadtxt("../circuit.dat", usecols=(0,1,2), unpack=True)
t = circ_ts*inputfile.dT*1e9 #ns
I = circ_deltaQ*convFactor2  #Amps
U = circ_U*convFactor3       #Volts

#Load heatspike data
print "Loading arcbounds_original.dat..."
(ts, emit, sigma, incident) = np.loadtxt("../arcbounds_original.dat", usecols=(0,10,11,12), unpack=True)
useIdxs_heatspike = np.arange(0,len(ts),5) ## !!! ASSUMPTION: 5 timesteps per heatspike !!!
ts = ts[useIdxs_heatspike]
emit = emit[useIdxs_heatspike]
sigma = sigma[useIdxs_heatspike]
incident = incident[useIdxs_heatspike]
sigma_idx = np.asarray(sigma,dtype='int')

t_hs = ts*inputfile.dT*1e9 #ns    
emit *= convFactor1/5.0 #particles/sec
sigma *= inputfile.dZ #cm
incident *= convFactor1/5.0 #particles/sec
area = np.pi*(sigma*1e4)**2 #um^2

# print t_hs
# print emit
# print sigma
# print incident
#exit(1)

#Load ez(t) data
fieldThreshold=-2000 #MV/m
rFactor     = 1e4*inputfile.Ldb #output units -> um
fieldFactor = 510.998928e3/2.99792458e10**2 * inputfile.Ldb * inputfile.O_pe**2 * 1e-4 #output units -> MV/m
#The data
field_ts = []
field_time = [] #ns
field_radius = [] #um
field_Ezave_rField = [] #MV/m
field_Ezave_rHS_time = [] #ns
field_Ezave_rHS = [] #MV/m
#Get list of output files and timestamp
timeIdx = open("../out/timeIndex.dat")
print "Reading Ez files..."
TIMEIDX_lineCounter = 0
idx_hspike = 0
timeIdx.readline() #Skip first line
for line in timeIdx:
    ls = line.split()
    stepNum=ls[0]
    timestamp = float(ls[1]) #ns

    if TIMEIDX_lineCounter % 10 == 0:
        print "stepNum =", stepNum
    TIMEIDX_lineCounter += 1
    
    fname = "../out/Ez" + stepNum + ".dat"
    fieldFile = open(fname, 'r')
    rList = []
    fieldList = []
    for line in fieldFile:
        l = line.split()
        if len(l) == 0:
            continue
        z = float(l[1])
        if (z > 0.0):
            #Want field at cathode only
            continue
        rList.append(float(l[0])*rFactor)
        fieldList.append(float(l[2])*fieldFactor)
    fieldFile.close()

    #Find point where field rises past 2000 MV/m, also calculate average field under sheath
    for i in xrange(idx_hspike,len(ts)):#start search at prev. known location, else it takes forever...
        if int(ts[i]) == int(stepNum):
            idx_hspike = i
            break
    fieldSum = fieldList[0]
    print sigma_idx[idx_hspike]
    for i in xrange(1,len(rList)):
        fieldSum += fieldList[i]
        if fieldList[i-1] < fieldThreshold and fieldList[i]>fieldThreshold:
            dFdr = (fieldList[i]-fieldList[i-1])/(rList[i]-rList[i-1])
            rField = rList[i-1]+(fieldThreshold-fieldList[i-1])/dFdr
            field_ts.append(stepNum)
            field_time.append(timestamp)
            field_radius.append(rField)
            
            field_Ezave_rField.append(fieldSum/(i+1))
        if i == int(sigma_idx[idx_hspike]):
            field_Ezave_rHS.append(fieldSum/(i+1))
            field_Ezave_rHS_time.append(t_hs[idx_hspike])
            print t_hs[idx_hspike], fieldSum/(i+1)
        
timeIdx.close()
field_ts = np.asarray(map(int,field_ts))
field_time = np.asarray(field_time)
field_radius = np.asarray(field_radius)
field_area = np.pi*field_radius**2

# print ts
# print field_ts
tsRatio = (field_ts[1]-field_ts[0])/(ts[1]-ts[0])
print "tsRatio =", tsRatio
#tsRatio = int(tsRatio)

def smooth(data,steps=1000, endFix=True):
    "Boxcar window smoother"
    print "Smoothing ("+str(steps) + "steps)..",
    if smooth.kernels==None:
        smooth.kernels = {}
    if not (steps in smooth.kernels):
        print "Creating kernel with steps =", \
            steps#, "=", steps*inputfile.dT*1e9, "ns =", steps*inputfile.dT*1e12, "ps", #steplength != sim timestep
        smooth.kernels[steps] = np.ones(steps)/float(steps)
    if len(data) < steps:
        print "Error: len(data) < steps!"
        exit(1)
    ret = np.convolve(data,smooth.kernels[steps],'same')
    if endFix:
        overhang1 = int(steps)/2
        overhang2 = steps-overhang1-1
        for i in xrange(overhang1):
            ret[i]    *= float(steps) / float( steps - overhang1 + i )
        for i in xrange(overhang2):
            ret[-i-1] *= float(steps) / float( steps - overhang2 + i )
    print "done."
    return ret
smooth.kernels=None
smoothsteps_hspike = 1000
smoothsteps_field = 10#int(smoothsteps_hspike/tsRatio)

### PLOTS ###
plt.figure(1)
plt.plot(sigma*1e4,I[useIdxs_heatspike], '.')
plt.xlabel("Radius [um]")
plt.ylabel("Current [A]")

plt.figure(2)
plt.plot(area,I[useIdxs_heatspike],'.',label="HS")
plt.xlabel("Area [um$^2$]")
plt.ylabel("Current [A]")

print "fitting I = a*AREA + b:"
cov = np.cov(area,I[useIdxs_heatspike])
fit_a = cov[0,1]/cov[0,0]
fit_b = np.mean(I[useIdxs_heatspike])-fit_a*np.mean(area)
fit_X = np.linspace(0.0,max(area))
fit_Y = fit_a*fit_X+fit_b
plt.plot(fit_X,fit_Y, lw=4,zorder=-1,color="red",label="Fit 1")
print "fit_a =", fit_a, "A/um^2, fit_b =", fit_b, "A"
fit_A2 = np.mean(area*I[useIdxs_heatspike])/np.mean(area**2)
print "fit_A2 =", fit_A2, "A/um^2"
fit_Y = fit_A2*fit_X
plt.plot(fit_X,fit_Y, lw=4,zorder=-1,color="green",ls="--", label="Fit 2")
plt.legend(loc=0)

plt.figure(3)
plt.plot(t_hs,sigma*1e4,label="HS")
plt.plot(t_hs,smooth(sigma*1e4,smoothsteps_hspike),label="HS smoothed",lw=4)
plt.plot(field_time, field_radius,label="Field, thresh="+str(fieldThreshold)+" MV/m")
plt.plot(field_time, smooth(field_radius,smoothsteps_field),label="Field, thresh="+str(fieldThreshold)+" MV/m (smoothed)", lw=4)
plt.xlabel("Time [ns]")
plt.ylabel("Radius [um]")
plt.legend(loc=0)

plt.figure(4)
plt.plot(t_hs,area,label="HS")
plt.plot(t_hs,smooth(area,smoothsteps_hspike),label="HS smoothed",lw=4)
plt.plot(field_time, field_area,label="Field, thresh="+str(fieldThreshold)+" MV/m")
plt.plot(field_time, smooth(field_area,smoothsteps_field),label="Field, thresh="+str(fieldThreshold)+" MV/m (smoothed)",lw=4)
plt.xlabel("Time [ns]")
plt.ylabel("Area [um$^2$]")
plt.legend(loc=0)

plt.figure(5)
plt.plot(field_area,I[field_ts], '+')
print "fitting I = a*AREA + b:"
cov = np.cov(field_area,I[field_ts])
fit_a = cov[0,1]/cov[0,0]
fit_b = np.mean(I[field_ts])-fit_a*np.mean(field_area)
fit_X = np.linspace(0.0,max(field_area))
fit_Y = fit_a*fit_X+fit_b
plt.plot(fit_X,fit_Y, lw=4,zorder=-1,color="red")
print "fit_a =", fit_a, "A/um^2, fit_b =", fit_b, "A"
fit_A2 = np.mean(field_area*I[field_ts])/np.mean(field_area**2)
print "fit_A2 =", fit_A2, "A/um^2"
fit_Y = fit_A2*fit_X
plt.plot(fit_X,fit_Y, lw=4,zorder=-1,color="green",ls="--")
plt.xlabel("Field emission area [um$^2$]")
plt.ylabel("Current [A]")

plt.figure(6)
plt.plot(t_hs,I[useIdxs_heatspike]/area,label="HS")
plt.plot(field_time,I[field_ts]/field_area,label="field")
plt.xlabel("Time [ns]")
plt.ylabel("Current density [A/um$^2$]")
plt.legend(loc=0)

plt.figure(7)
plt.plot(t_hs,smooth(I[useIdxs_heatspike]/area,smoothsteps_hspike),label="HS")
plt.plot(field_time,smooth(I[field_ts]/field_area,smoothsteps_field),label="field")
plt.xlabel("Time [ns]")
plt.ylabel("Current density [A/um$^2$]")
plt.legend(loc=0)

plt.figure(8)
plt.plot(field_time, field_Ezave_rField, label="r from field")
plt.plot(field_Ezave_rHS_time, field_Ezave_rHS, label="r from HS")
plt.xlabel("Time [ns]")
plt.ylabel("Average field under sheath [MV/m]")
plt.legend(loc=0)

plt.show()
