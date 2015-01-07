#!/usr/bin/env python
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
# plot_compare.py:
# Plot currents, particle numbers and more. Compare several runs in one plot.
#

import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rcParams,rc
rcParams.update({'text.usetex': True})
PLOTOPTS = 1
if PLOTOPTS == 1: #Refrun figure for paper
    DPI = 500
    rcParams.update({'savefig.dpi':DPI})
    rc('font',**{'family':'serif','serif':['Times'],'size':8})
    textwidth = 449.40363*1.0/72.27 #inches, for 2Dpic paper
    FIGSIZE = (0.5*0.95*textwidth, 0.5*0.98*textwidth/1.61803398875)
    altLS = "--"
elif PLOTOPTS == 2 or PLOTOPTS == 3: #surface model figure or circuit params figure for paper
    DPI = 500
    rcParams.update({'savefig.dpi':DPI})
    rc('font',**{'family':'serif','serif':['Times'],'size':8})
    textwidth = 449.40363*1.0/72.27 #inches, for 2Dpic paper
    FIGSIZE = (0.5*0.95*textwidth, 0.5*0.98*textwidth/1.61803398875)
    altLS = ":"
else:
    FIGSIZE = None
    DPI=None
    altLS = "--"

print "PLOTOPTS =", PLOTOPTS

import sys

from calcScaling import InputFile


def Usage():
    print "plot_compareIV.py <IV|hspike|nparticles> folderName--title folderName--title ..."
    exit(1)
if len(sys.argv) <= 1:
    Usage()

mode = sys.argv[1]
assert mode == "IV" or mode == "hspike" or mode=="nparticles"

def smooth(data,steps, endFix=True):
    "Boxcar window smoother"
    if smooth.kernel==None or len(smooth.kernel) != steps:
        print "Creating kernel with steps =", steps
        smooth.kernel = np.ones(steps)/float(steps)
    if len(data) < steps:
        print "Error: len(data) < steps!"
        exit(1)
    ret = np.convolve(data,smooth.kernel,'same')
    if endFix:
        overhang1 = int(steps)/2
        overhang2 = steps-overhang1-1
        for i in xrange(overhang1):
            ret[i]    *= float(steps) / float( steps - overhang1 + i )
        for i in xrange(overhang2):
            ret[-i-1] *= float(steps) / float( steps - overhang2 + i )
    return ret
smooth.kernel=None

def getData_IV(folderName):
    print "********* Loading IV data from '" + folderName + "' ***********"
    inputfile = InputFile(folderName + "/input.txt")
    inputfile.calcBasicParams()
    
    convFactor1 = inputfile.N_sp/inputfile.dT #superpartices/timestep -> particles/sec
    convFactor2 = 1.60217657e-19*inputfile.N_sp/inputfile.dT # superparticles/timestep -> Amps
    convFactor3 = (inputfile.dz/inputfile.Omega_pe)**2*inputfile.T_ref #dimless potential -> Volts
    
    #Load circuit data
    if inputfile.CircuitName == "FixedVoltage_resistorCapacitor":
        (circ_ts, circ_deltaQ, circ_U, circuitCurrent) = \
            np.loadtxt(folderName + "/circuit.dat", usecols=(0,1,2,3), unpack=True)
        
        t = circ_ts*inputfile.dT*1e9         #ns
        I = circ_deltaQ*convFactor2          #Amps
        U = circ_U*convFactor3               #Volts
        I_circ = circuitCurrent*convFactor2  #Amps
        
        print "done."
        return (t,I,U,I_circ)
    
    else:
        (circ_ts, circ_deltaQ, circ_U) = np.loadtxt(folderName + "/circuit.dat", usecols=(0,1,2), unpack=True)
        circuitCurrent=np.zeros_like(circ_ts)

        t = circ_ts*inputfile.dT*1e9         #ns
        I = circ_deltaQ*convFactor2          #Amps
        U = circ_U*convFactor3               #Volts
    
        print "done."
        return (t,I,U)

def getDATA_hspike(folderName):
    print "********* Loading hspike data from '" + folderName + "' ***********"
    inputfile = InputFile(folderName + "/input.txt")
    inputfile.calcBasicParams()

    convFactor1 = inputfile.N_sp/inputfile.dT #superpartices/timestep -> particles/sec
    convFactor2 = 1.60217657e-19*inputfile.N_sp/inputfile.dT # superparticles/timestep -> Amps
    convFactor3 = (inputfile.dz/inputfile.Omega_pe)**2*inputfile.T_ref #dimless potential -> Volts
    
    (ts, emit, sigma, incident) = np.loadtxt(folderName + "/arcbounds_original.dat", usecols=(0,10,11,12), unpack=True)
    useIdxs = np.arange(0,len(ts),5)
    ts = ts[useIdxs]
    emit = emit[useIdxs]
    sigma = sigma[useIdxs]
    incident = incident[useIdxs]
    
    t = ts*inputfile.dT*1e9 #ns    
    emit *= convFactor1/5.0 #particles/sec
    sigma *= inputfile.dZ #cm
    incident *= convFactor1/5.0 #particles/sec
    
    return (t,emit,sigma,incident)

def getDATA_nParticles(folderName):
    print "********* Loading nParticles data from '" + folderName + "' ***********"
    inputfile = InputFile(folderName + "/input.txt")
    inputfile.calcBasicParams()

    (t, e,i,n) = np.loadtxt(folderName + "/mainStats.dat", usecols=(1,3,4,5), unpack=True)
    
    e*=inputfile.N_sp
    i*=inputfile.N_sp
    n*=inputfile.N_sp
    
    return (t,e,i,n)

FOLDERS  = []
DATA  = []
TITLES   = []
V0_all = []
for i in xrange(2,len(sys.argv)):
    foo = sys.argv[i]
    (folder,title) = foo.split("--")
    
    FOLDERS.append(folder)
    TITLES.append(title)
    
    print "Now plotting title='"+title+"'"
    
    if mode == "IV":
        DATA.append(getData_IV(folder))

        plt.figure(1,figsize=FIGSIZE,dpi=DPI) #Current
        currentLine, = plt.plot(DATA[-1][0],DATA[-1][1], label=title)
        if len(DATA[-1]) == 4: # has circuitCurrent
            plt.plot(DATA[-1][0],DATA[-1][3], ls=altLS,color=currentLine.get_color())

        plt.figure(2,figsize=FIGSIZE,dpi=DPI) #Voltage
        plt.plot(DATA[-1][0],DATA[-1][2], label=title)
        
        if PLOTOPTS != 3:
            #Per-simulation plots
            plt.figure(3+(i-2),figsize=FIGSIZE,dpi=DPI) #IV / single run
            plt.title(title)
            line_Igap, = plt.plot(DATA[-1][0],DATA[-1][1], 'b-')
            line_Icirc = None
            if len(DATA[-1]) == 4: # has circuitCurrent
                line_Icirc, = plt.plot(DATA[-1][0],DATA[-1][3], 'g-')
            #plt.plot(DATA[-1][0],DATA[-1][1], 'b-')
            ax1 = plt.gca()
            ax1.set_ylabel("Current [A]", color='b')
            ax1.set_xlabel("Time [ns]")
            for tl in ax1.get_yticklabels():
                tl.set_color('b')

            ax2 = ax1.twinx()
            if PLOTOPTS == 1:
                line_Vgap, = ax2.plot(DATA[-1][0],DATA[-1][2]*1e-3, 'r-')
                ax2.set_ylabel("Voltage [kV]", color='r')
            else:
                line_Vgap, = ax2.plot(DATA[-1][0],DATA[-1][2], 'r-')
                ax2.set_ylabel("Voltage [V]", color='r')
            for tl in ax2.get_yticklabels():
                tl.set_color('r')
            ax2.yaxis.get_offset_text().set_color('r')

            if line_Igap != None:
                ax2.legend((line_Icirc, line_Igap, line_Vgap),(r"$I_\mathrm{circ}$",r"$I_\mathrm{gap}$",r"$V_\mathrm{gap}$"), loc=0,frameon=False)

            if PLOTOPTS == 1:
                #ax2.get_yaxis().get_major_formatter().set_powerlimits((0, 0))
                #plt.subplots_adjust(right=0.84, left=0.1, top=0.92, bottom=0.17) #1e3 V
                plt.subplots_adjust(right=0.84, left=0.1, top=0.97, bottom=0.17) #kV
                plt.xlim(0,2.1)


        #find t90 and t10 for risetime calc
        V0  = DATA[-1][2][0]
        V0_all.append(V0)
        V90 = DATA[-1][2][0]*0.9
        V50 = DATA[-1][2][0]*0.5
        V10 = DATA[-1][2][0]*0.1
        t90 = None
        t50 = None
        t10 = None
        for i in xrange(0, len(DATA[-1][0])):
            if t90 == None and DATA[-1][2][i] <= V90:
                t90 = DATA[-1][0][i]
            elif t50 == None and DATA[-1][2][i] <= V50:
                t50 = DATA[-1][0][i]
            elif t10 == None and DATA[-1][2][i] <= V10:
                t10 = DATA[-1][0][i]
                break
        print "V0  =", V0,  "[V]"
        print "V90 =", V90, "[V], t90 =", t90, "[ns]"
        print "V10 =", V10, "[V], t10 =", t10, "[ns]"
        
        if PLOTOPTS == 3:
            continue
        
        if t90 != None:
            plt.axvline(t90, ls="--", color="k")
        if t50 != None:
            print "t50-t90 =", t50-t90, "[ns]"
            plt.axvline(t50, ls="--", color="k")
            if PLOTOPTS == 1:
                plt.annotate(
                    '', xy=(t90, V50*1e-3), xycoords = 'data',
                    xytext = (t50, V50*1e-3), textcoords = 'data',
                    arrowprops = {'arrowstyle':'<->'})
                plt.annotate(
                    ("$%.2f" % (t50-t90) ) + r'$',
                    xy=(t90 + (t50-t90)/2.0, 1e-3*(V50+50)), xycoords = 'data',
                    horizontalalignment='center', verticalalignment='bottom')
            else:
                plt.annotate(
                    '', xy=(t90, V50), xycoords = 'data',
                    xytext = (t50, V50), textcoords = 'data',
                    arrowprops = {'arrowstyle':'<->'})
                plt.annotate(
                    ("$%.2f" % (t50-t90) ) + r'$ [ns]',
                    xy=(t90 + (t50-t90)/2.0, V50), xycoords = 'data',
                    horizontalalignment='center', verticalalignment='bottom')
        if t10 != None:
            print "t10-t90 =", t10-t90, "[ns]"
            plt.axvline(t10, ls="--", color="k")

            if PLOTOPTS == 1:
                plt.annotate(
                    '', xy=(t10, V90*1e-3), xycoords = 'data',
                    xytext = (t90, V90*1e-3), textcoords = 'data',
                    arrowprops = {'arrowstyle':'<->'})
                plt.annotate(
                    ("$%.2f" % (t10-t90) ) + r'$',
                    xy=(t90 + (t50-t90)/2.0, 1e-3*(V90+50)), xycoords = 'data',
                    horizontalalignment='center', verticalalignment='bottom')
            else:
                plt.annotate(
                    ("$%.2f" % (t10-t90) ) + r'$ [ns]',
                    xy=(t90 + (t10-t90)/2.0, V90), xycoords = 'data',
                    horizontalalignment='center', verticalalignment='bottom')
                plt.annotate(
                    '', xy=(t10, V90), xycoords = 'data',
                    xytext = (t90, V90), textcoords = 'data',
                    arrowprops = {'arrowstyle':'<->'})
        if PLOTOPTS == 1:
            #For B/W compatability
            plt.annotate(
                r"$V_\textrm{gap}$", xy = (0.5,1740*1e-3), xycoords='data',
                color='red', xytext=(0, -20), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),
                horizontalalignment='center', verticalalignment='top'
                )
            plt.annotate(
                r"$I_\textrm{gap}$", xy = (1.25,-300*1e-3), xycoords='data',
                color='blue', xytext=(-20, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),
                horizontalalignment='center', verticalalignment='bottom'
                )
            plt.annotate(
                r"$I_\textrm{circ}$", xy = (1.5,-350*1e-3), xycoords='data',
                color='green', xytext=(0, 15), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),
                horizontalalignment='center', verticalalignment='bottom'
                )
            
        print
        print

    elif mode == "hspike":
        DATA.append(getDATA_hspike(folder))
        
        plt.figure(1) #radius
        plt.plot(DATA[-1][0], DATA[-1][2]*1e4, label=title)
        
        plt.figure(2)
        plt.plot(DATA[-1][0], smooth(DATA[-1][2]*1e4, 1000), label=title)

        plt.figure(3) #emitted
        plt.plot(DATA[-1][0], DATA[-1][1], label=title)

        plt.figure(4) #incident
        plt.plot(DATA[-1][0], DATA[-1][3], label=title)
        
        plt.figure(5)
        plt.plot(DATA[-1][0], DATA[-1][3]/(np.pi*DATA[-1][2]**2), label=title)

        plt.figure(6)
        plt.plot(DATA[-1][0], smooth(DATA[-1][3]/(np.pi*DATA[-1][2]**2),1000), label=title)        
    
    elif mode == "nparticles":
        DATA.append(getDATA_nParticles(folder))
        
        plt.figure(1)
        plt.plot(DATA[-1][0], DATA[-1][1],label=title)
        
        plt.figure(2)
        plt.plot(DATA[-1][0], DATA[-1][2],label=title)
        
        plt.figure(3)
        plt.plot(DATA[-1][0], DATA[-1][3],label=title)
    
    print
    print
    print
    
if mode == "IV":
    plt.figure(1)
    plt.xlabel("Time [ns]")
    plt.ylabel("Current [A]")
    if PLOTOPTS == 1:
        plt.subplots_adjust(right=0.88, left=0.15, top=0.92, bottom=0.17)
        plt.xlim(0,2.1)
        plt.legend(loc=2,frameon=False)
    elif PLOTOPTS == 2:
        plt.subplots_adjust(right=0.97, left=0.11, top=0.96, bottom=0.16)
        plt.legend(loc=2,frameon=False,fontsize=8)
        plt.xlim(0,2.5)
    elif PLOTOPTS == 3:
        plt.subplots_adjust(right=0.97, left=0.12, top=0.96, bottom=0.16)
        plt.legend(loc=2,frameon=False,fontsize=8)
        #plt.xlim(0,2.5)
    else:
        plt.legend(loc=2,frameon=False)

    plt.figure(2)
    plt.xlabel("Time [ns]")
    plt.ylabel("Voltage [V]")
    if PLOTOPTS == 1:
        plt.subplots_adjust(right=0.88, left=0.15, top=0.92, bottom=0.17)
        plt.xlim(0,2.1)
        plt.legend(loc=3,frameon=False)
    elif PLOTOPTS == 2:
        plt.subplots_adjust(right=0.97, left=0.17, top=0.96, bottom=0.16)
        plt.xlim(0,2.5)
        plt.ylim(-750,2000)
        plt.legend(loc=3,frameon=False,fontsize=8)
    elif PLOTOPTS == 3:
        plt.subplots_adjust(right=0.97, left=0.17, top=0.96, bottom=0.16)
        #plt.xlim(0,2.5)
        #plt.ylim(-750,2000)
        plt.legend(loc=3,frameon=False,fontsize=8)
    else:
        plt.legend(loc=3,frameon=False)
    V0 = V0_all[0]
    allV0EqualFlag = True
    for v0 in V0_all:
        if v0 != V0:
            allV0EqualFlag == False
            break
    if allV0EqualFlag:
        if PLOTOPTS == 2:
            plt.axhline(0.9*V0,color='k',ls='--',xmin=0.4)
            plt.axhline(0.5*V0,color='k',ls='--',xmin=0.5)
            plt.axhline(0.1*V0,color='k',ls='--',xmin=0.6)
        if PLOTOPTS == 3:
            plt.axhline(0.9*V0,color='k',ls='--',xmin=0.0)
            plt.axhline(0.5*V0,color='k',ls='--',xmin=0.5)
            plt.axhline(0.1*V0,color='k',ls='--',xmin=0.5)
        else:
            plt.axhline(0.9*V0,color='k',ls='--')
            plt.axhline(0.5*V0,color='k',ls='--')
            plt.axhline(0.1*V0,color='k',ls='--')


elif mode == "hspike":
    plt.figure(1)
    plt.xlabel("Time [ns]")
    plt.ylabel("Radius [um]")
    plt.legend(loc=0)

    plt.figure(2)
    plt.xlabel("Time [ns]")
    plt.ylabel("Radius [um]")
    plt.legend(loc=0)
    #plt.xlim(0,0.4) ##

    plt.figure(3)
    plt.xlabel("Time [ns]")
    plt.ylabel("Particles/sec")
    plt.legend(loc=0)

    plt.figure(4)
    plt.xlabel("Time [ns]")
    plt.ylabel("Particles/sec")
    plt.legend(loc=0)

    plt.figure(5)
    plt.xlabel("Time [ns]")
    plt.ylabel("Particles/cm$^2$/sec")
    plt.legend(loc=0)

    plt.figure(6)
    plt.xlabel("Time [ns]")
    plt.ylabel("Particles/cm$^2$/sec")
    plt.legend(loc=0)
elif mode == "nparticles":
    plt.figure(1)
    plt.xlabel("Time [ns]")
    plt.ylabel("Number of electrons")
    plt.legend(loc=0)
    (ymin,ymax) = plt.ylim()
    the_ymax = ymax
    #plt.xlim(0,0.4) ##

    plt.figure(2)
    plt.xlabel("Time [ns]")
    plt.ylabel("Number of ions")
    plt.legend(loc=0)
    (ymin,ymax) = plt.ylim()
    if the_ymax < ymax:
        the_ymax = ymax
    #plt.xlim(0,0.4) ##

    plt.figure(3)
    plt.xlabel("Time [ns]")
    plt.ylabel("Number of neutrals")
    plt.legend(loc=0)
    (ymin,ymax) = plt.ylim()
    if the_ymax < ymax:
        the_ymax = ymax
    #plt.xlim(0,0.4) ##

    plt.figure(1)
    plt.ylim(ymax=the_ymax)
    plt.figure(2)
    plt.ylim(ymax=the_ymax)
    plt.figure(3)
    plt.ylim(ymax=the_ymax)
    
plt.show()
    
