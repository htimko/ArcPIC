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
# particleHisto.py:
# Calculates the particle density from the r{e|Cu|Cup}...dat files
# using a histogram, and writes the results to a new file.
# This can then be used to plot 2D particle density maps.
#

import sys, os
import numpy as np

if not len(sys.argv) == 11:
    print "Usage: ./particleHisto.py <inFile> <outfile> zBins rBins zMax rMax dZ dz Rsp getLog"
    exit(1)

inFileName  = sys.argv[1]
outFileName = sys.argv[2]
zBins = int(sys.argv[3])
rBins = int(sys.argv[4])
zMax = float(sys.argv[5])    # system size [dz]
rMax = float(sys.argv[6])
dZ = float(sys.argv[7]) # [um/grid]
dz = float(sys.argv[8]) # [grid/Ldb] (typically 0.5)
N_sp  = float(sys.argv[9])   # num real particles / superparticle
getLog = False
if sys.argv[10] == "y":
    getLog = True
elif sys.argv[10] == "n":
    getLog = False
else:
    print "Error: getLog must be 'y' or 'n'!"
    exit(1)

dz_bin = zMax/float(zBins)
dr_bin = rMax/float(rBins)

print "Using binsize (dz_bin,dr_bin) = (%g,%g)" % (dz_bin,dr_bin);

histo = np.zeros((rBins,zBins),dtype=np.float64)

if os.path.isfile(outFileName):
    print "Error: Can't create output file, '" + outFileName + "' already exists!"
    exit(1)

#Fill histogram
inFile = open(inFileName, 'r');
for line in inFile:
    l = line.split()
    z  = float(l[0])/dz; r  = float(l[1])/dz;
    zi = int(z/dz_bin);   ri = int(r/dr_bin)

    if zi >= zBins or ri >= rBins:
        print "Error: (zi,ri) = (%i,%i) out of bounds (%i,%i)" % (zi, ri, zBins, rBins)
        exit(1)
    
    histo[ri, zi] += 1;
inFile.close()

#Normalize to num. physical particles / cm^2
histo[:,:] *= N_sp
for ri in xrange(0,rBins):
    histo[ri,:] /= np.pi*(2*ri+1)*(dZ*1e-4)**3

#output to file
outFile = open(outFileName, 'w');
outFile.write("! nx %i ny %i xmin 0.0 xmax %g ymin 0.0 ymax %g\n" % (zBins, rBins, zMax*dZ, rMax*dZ))
for ri in xrange(0,rBins):
    for zi in xrange(0,zBins):
        if getLog:
            v = histo[ri,zi]
            if v > 0.0:
                v = np.log10(v)
            else:
                v = 0.0
            outFile.write("%g " % (v,) );
        else:
            outFile.write("%g " % (histo[ri,zi],) );
    outFile.write("\n");
outFile.close()
