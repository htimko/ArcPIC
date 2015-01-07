#!/usr/bin/env python

import matplotlib.pyplot as plt
import os
import re
import math as m
import numpy as np

import matplotlib as mpl
mpl.rcParams['text.usetex'] = True


#Physics settings
Te = 4.0 #eV
ne = 10e10; #cm^-3

c = 3e8
me = 511e3 / c**2 #eV/c^2

vte = m.sqrt(Te/me) #Thermal velocity [m/s]
vk  = m.sqrt(2)*vte   #Kinetic velocity (most probable value) [m/s]

#nu_0 = 6.41306e-32 * 177.715 * (Te*8.85418e-12)**(3.0/2.0) / \
#     ( m.sqrt(ne)* (1.60217646e-19)**3.0 * 13.0 )
nu_0 = 0.0000502601*ne/m.sqrt(Te**3)
omega_pe = 56414.6*m.sqrt(ne);
steps_nu_0 = omega_pe/(0.2*nu_0)
print "nu_0 =", nu_0, " [s^-1], omega_pe =", omega_pe, "[s^-1], nu_0^-1/DeltaT =", steps_nu_0, "[steps]"

def Maxwell(v):
    "Maxwell distribution"
    return (me/(2*m.pi*Te))**(3.0/2.0) * 4*m.pi * v**2 * np.exp(-me*v**2/(2*Te))
def Gaus(v,mu,sigma):
    "Gaussian distribution"
    return np.exp(-0.5*((v-mu)/sigma)**2)/(sigma*m.sqrt(2*m.pi))

#Read the energy file
efile = open("energy.out",'r');
eTot = []; eMean = [];
for line in efile:
    (i, et,em) = line.split()
    eTot.append(et);
    eMean.append(em);
efile.close()

plt.figure(1);
plt.plot(eMean)
#plt.show()

#Get the velocity files
fileList = os.listdir(os.getcwd());
veloFiles = []
for f in fileList:
    if f.startswith("vdist_") and f.endswith(".out"):
        veloFiles.append(f);
veloFiles = sorted(veloFiles)

#Loop over velocity files
mybins_v = np.linspace(0,vk*4,50)
mybins_vxyz = np.linspace(-vk*4/m.sqrt(3),vk*4/m.sqrt(3),50)
for i in xrange(len(veloFiles)):
    stepNum = int(re.match("vdist_(\d+).out", veloFiles[i]).group(1))
    
    #Read the file
    vfile = open(veloFiles[i],'r')
    vx = []; vy = []; vz = []; v = [];
    for line in vfile:
        (x,y,z) = line.split(", ")
        x = float(x); y = float(y); z = float(z)
        vx.append(x)
        vy.append(y)
        vz.append(z)
        
        v.append(m.sqrt(x**2 + y**2 + z**2))

        #print v[-1] / 1.18701e+06;
    
    #Plot histograms
    plt.figure(2);
    (pdf, bins, patches) = plt.hist(v, bins=mybins_v, normed=True, alpha=0.3, label = r'$\nu_0~t = %f$' % (stepNum/float(steps_nu_0),), histtype="stepfilled" )
    #print pdf, bins
    #print np.sum(pdf * np.diff(bins))

    plt.figure(3)
    (pdf, bins, patches) = plt.hist(vx, bins=mybins_vxyz, normed=True, alpha=0.3, label = r'$\nu_0~t = %f$' % (stepNum/float(steps_nu_0),), histtype="stepfilled" )
    plt.figure(4)
    (pdf, bins, patches) = plt.hist(vy, bins=mybins_vxyz, normed=True, alpha=0.3, label = r'$\nu_0~t = %f$' % (stepNum/float(steps_nu_0),), histtype="stepfilled" )
    plt.figure(5)
    (pdf, bins, patches) = plt.hist(vz, bins=mybins_vxyz, normed=True, alpha=0.3, label = r'$\nu_0~t = %f$' % (stepNum/float(steps_nu_0),), histtype="stepfilled" )

plt.figure(2)
vaxis = np.linspace(0,vk*4)
plt.plot(vaxis, Maxwell(vaxis), 'r', label="Maxwell");
plt.legend()

vxyzaxis = np.linspace(-vk*4/m.sqrt(3),vk*4/m.sqrt(3))
plt.figure(3)
plt.plot(vxyzaxis,Gaus(vxyzaxis,0.0,m.sqrt(Te/me)), label="Maxwell")
plt.legend()
plt.figure(4)
plt.plot(vxyzaxis,Gaus(vxyzaxis,0.0,m.sqrt(Te/me)), label="Maxwell")
plt.legend()
plt.figure(5)
plt.plot(vxyzaxis,Gaus(vxyzaxis,0.0,m.sqrt(Te/me)), label="Maxwell")
plt.legend()

plt.show()
        
