#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# plotPicFieldTest.py
#
# Plots the PIC and direct field as well as their difference
# in r=z=0 as function of charge position
#
# Kyrre Sjøbæk, 2013
#

import sys, os, shutil
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
#from matplotlib.ticker import LogFormatter

#from matplotlib import gridspec

from matplotlib import rcParams
rcParams.update({'text.usetex': True})

print "Usage: {g|ng} {mirrorOrders=y|n} {ZR=z|r} (maxZidx|None) (hlineValue)"
print
if (len(sys.argv) < 4) != (len(sys.argv) > 6):
    exit(1)


#Get the scaling setup
MODLOAD_parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,MODLOAD_parentdir) 
from calcScaling import InputFile
inputfile = InputFile("../input.txt")
inputfile.calcBasicParams()
print

rFactor   = 1e4*inputfile.dZ #output units -> um

idxList = []
izList  = []
irList  = []
rList   = []
zList   = []
fieldPICList = []
fieldDirectList = []

fieldDirect_orderList = []


useGrid = sys.argv[1]
if useGrid == 'g':
    useGrid = True
elif useGrid == 'ng':
    useGrid = False
else:
    print "argument 1 (useGrid) should be 'g' or 'ng'"
    exit(1)

mirrorOrders = sys.argv[2]
if mirrorOrders == 'y':
    mirrorOrders = True
elif mirrorOrders == 'n':
    mirrorOrders = False
else:
    print "argument 2 (mirrorOrders) should be 'y' or 'n'"
    exit(1)

ZR = sys.argv[3]
if not (ZR == 'z' or ZR == 'r'):
    print "Argument 3 (ZR) should be 'z' or 'r', got '" + ZR + "'"
    exit(1)

maxZidx = None
if len(sys.argv) == 5 or len(sys.argv) == 6:
    if sys.argv[4] != "None":
        maxZidx = int(sys.argv[4])
        print "Using maxZidx =", maxZidx
        print

hlineValue = None
if len(sys.argv) == 6:
    hlineValue = float(sys.argv[5])
    print "Using hlineValue =", hlineValue
    print

picFieldTestFile = open("picFieldTest.dat", 'r')
picFieldTestFile.readline() #skip header

header2 = picFieldTestFile.readline()
h2 = header2.split()
assert len(h2) == 7
nzMax = int(h2[1].split("=")[1])
nrMax = int(h2[2].split("=")[1])
pointsPerGrid = int(h2[3].split("=")[1])
mirror_order = int(h2[4].split("=")[1])
zWitness = float(h2[5].split("=")[1])
rWitness = float(h2[6].split("=")[1])


for line in picFieldTestFile:
    l = line.split()
    if len(l) == 0:
        continue
    
    if maxZidx != None and int(l[1]) >= maxZidx:
        print "stop"
        break;
    
    idxList.append(           int(l[0]) )
    
    izList.append(            int(l[1]) )
    irList.append(            int(l[2]) )

    zList.append(           float(l[3]) )
    rList.append(           float(l[4]) )
    
    if ZR == 'z':
        fieldPICList.append(    float(l[5]) )
    elif ZR == 'r':
        fieldPICList.append(    float(l[6]) )
    
    fieldDirectList.append( float(l[7]) )
    
    fieldDirect_order = l[8:]
    fieldDirect_orderList.append(map(float,fieldDirect_order))

picFieldTestFile.close()

print "Max iz =", max(izList), "(pass this-1 as argument maxZidx if crash)"
print

idxList = np.asarray(idxList)
izList  = np.asarray(izList)
irList  = np.asarray(irList)
zList   = np.asarray(zList)
rList   = np.asarray(rList)
fieldPICList = np.asarray(fieldPICList)
fieldDirectList = np.asarray(fieldDirectList)

fieldDirect_orderList = np.asarray(fieldDirect_orderList)

zMat = np.reshape(zList, (-1,nrMax*pointsPerGrid))
rMat = np.reshape(rList, (-1,nrMax*pointsPerGrid))
fieldPICMat = np.reshape(fieldPICList, (-1, nrMax*pointsPerGrid))
fieldDirectMat = np.reshape(fieldDirectList, (-1, nrMax*pointsPerGrid))

def plotGrid():
    if not useGrid:
        return
    for ir in xrange(1, int(max(rList))+1 ):
        plt.axhline(ir, color='k', ls='--')
    for iz in xrange(1, int(max(zList))+1 ):
        plt.axvline(iz, color='k', ls='--')
    plt.plot(zWitness,rWitness, marker='*', ms=30)

fieldMax = max(np.nanmax(fieldPICMat), np.nanmax(fieldDirectMat))
fieldLevels_maxExp = int(np.log10(fieldMax))+1
fieldLevels_steps = fieldLevels_maxExp + 4 #starts at -3, all ints
fieldLevels_steps += fieldLevels_steps - 1  #add half-steps
fieldLevels = np.logspace(-3, fieldLevels_maxExp, fieldLevels_steps);
print "FieldLevels = "
print fieldLevels
print

fieldNegMax = max(np.nanmax(-fieldPICMat), np.nanmax(-fieldDirectMat))
fieldLevels_negMaxExp = int(np.log10(fieldNegMax))+1
fieldLevels_negSteps = fieldLevels_negMaxExp + 4 #starts at -3, all ints incl. 0
fieldLevels_negSteps += fieldLevels_negSteps - 1    #add half-steps
fieldLevels_neg = np.logspace(-3, fieldLevels_negMaxExp, fieldLevels_negSteps);
print "FieldLevels_neg = "
print fieldLevels_neg
print

plt.figure(1)
plt.contourf(zMat,rMat,fieldPICMat, 20);
plt.colorbar()
if hlineValue:
    plt.axhline(hlineValue, color='k')
plt.xlabel("z [dz=%f um]" % (inputfile.dZ*1e4,))
plt.ylabel("r [dz=%f um]" % (inputfile.dZ*1e4,))
plt.title("PIC field [V/m]")
plotGrid()

plt.figure(2)
CS1 = plt.contourf(zMat,rMat,fieldPICMat,
                   levels=fieldLevels, norm=LogNorm());
plt.colorbar(CS1)
CS2 = plt.contourf(zMat,rMat,-fieldPICMat,
                   levels=fieldLevels_neg, norm=LogNorm(), hatches='x');
plt.colorbar(CS2)
if hlineValue:
    plt.axhline(hlineValue, color='k')
plt.xlabel("z [dz=%f um]" % (inputfile.dZ*1e4,))
plt.ylabel("r [dz=%f um]" % (inputfile.dZ*1e4,))
plt.title("PIC field [V/m]")
plotGrid()

plt.figure(3)
plt.contourf(zMat,rMat,fieldDirectMat);
plt.colorbar()
if hlineValue:
    plt.axhline(hlineValue, color='k')
plt.xlabel("z [dz=%f um]" % (inputfile.dZ*1e4,))
plt.ylabel("r [dz=%f um]" % (inputfile.dZ*1e4,))
plt.title("Analytical field [V/m]")
plotGrid()

plt.figure(4)
CS1 = plt.contourf(zMat,rMat,fieldDirectMat,
                   levels=fieldLevels, norm=LogNorm());
plt.colorbar(CS1)
CS2 = plt.contourf(zMat,rMat,-fieldDirectMat,
                   levels=fieldLevels_neg, norm=LogNorm(), hatches='x');
plt.colorbar(CS2)
if hlineValue:
    plt.axhline(hlineValue, color='k')
plt.xlabel("z [dz=%f um]" % (inputfile.dZ*1e4,))
plt.ylabel("r [dz=%f um]" % (inputfile.dZ*1e4,))
plt.title("Analytical field [V/m]")
plotGrid()

deltaFieldMat = fieldPICMat - fieldDirectMat
plt.figure(5)
plt.contourf(zMat,rMat,deltaFieldMat);
plt.colorbar()
if hlineValue:
    plt.axhline(hlineValue, color='k')
plt.xlabel("z [dz=%f um]" % (inputfile.dZ*1e4,))
plt.ylabel("r [dz=%f um]" % (inputfile.dZ*1e4,))
plt.title("PIC - analytical field [V/m]")
plotGrid()

plt.figure(6)
dfmLevels_max = np.nanmax(abs(deltaFieldMat))
dfmLevels_maxExp = int(np.log10(dfmLevels_max))+1
dfmLevels_steps = dfmLevels_maxExp + 4 #starts at -3, all ints
dfmLevels_steps += dfmLevels_steps - 1 #add half-steps
dfmLevels = np.logspace(-3, dfmLevels_maxExp, dfmLevels_steps);
print "dfmLevels = "
print dfmLevels
print

plt.contourf(zMat,rMat,-deltaFieldMat,
             levels=dfmLevels,
             norm=LogNorm(),
             hatches="x");
plt.contourf(zMat,rMat,deltaFieldMat,
             levels=dfmLevels,
             norm=LogNorm());
if hlineValue:
    plt.axhline(hlineValue, color='k')
#plt.contourf(zMat,rMat,abs(deltaFieldMat), 20, norm=LogNorm());
plt.colorbar()
#plt.contour(zMat,rMat,deltaFieldMat, [0.0], norm=LogNorm());
plt.xlabel("z [dz=%f um]" % (inputfile.dZ*1e4,))
plt.ylabel("r [dz=%f um]" % (inputfile.dZ*1e4,))
plt.title("PIC - analytical field [V/m]\n(PIC underestimated region hatched)")
plotGrid()

relDeltaFieldMat = deltaFieldMat / abs(fieldDirectMat)
plt.figure(7)
#relDeltaFieldMat_levels = np.linspace(min(relDeltaFieldMat),max(relDeltaFieldMat), 10)
#relDeltaFieldMat_levels = np.arange(-1.0,1.1, 0.1)
relDeltaFieldMat_levels = (-1.0,-0.8,-0.6,-0.4,-0.2,-0.1,-0.05, 0.0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0)
print "relDeltaFieldMat_levels:"
print relDeltaFieldMat_levels
print
CS1 = plt.contourf(zMat,rMat,relDeltaFieldMat, relDeltaFieldMat_levels,extend='both');
plt.colorbar(CS1)
CS2 = plt.contour(zMat,rMat,relDeltaFieldMat, [-0.1,-0.05,0.0,0.05,0.1], colors='k');
plt.clabel(CS2,[-0.1,0.0,0.1])
if hlineValue:
    plt.axhline(hlineValue, color='k')
plt.xlabel("z [dz=%f um]" % (inputfile.dZ*1e4,))
plt.ylabel("r [dz=%f um]" % (inputfile.dZ*1e4,))
plt.title("(PIC - analytical field) / abs(analytical field) [V/m]\nBlack lines drawn at 0.0, $\pm$0.05 and $\pm$0.1")
plotGrid()

# plt.figure(8)
# plt.plot(zList,rList, 'b.')
# plt.xlabel("z [dz=%f um]" % (inputfile.dZ*1e4,))
# plt.ylabel("r [dz=%f um]" % (inputfile.dZ*1e4,))
# plotGrid()

def extractRslice(rWanted):
    indexList_wanted = []
    for i in xrange(len(rList)):
        if rList[i] == rWanted: #float stored and retrived without any operations in between
            indexList_wanted.append(i)
    zList_wanted = zList[indexList_wanted]
    fieldPICList_wanted = fieldPICList[indexList_wanted]
    fieldDirectList_wanted = fieldDirectList[indexList_wanted]
    return (zList_wanted, fieldPICList_wanted, fieldDirectList_wanted)

RValues = sorted(list(set(rList)))
print "RValues = ", RValues

plt.figure(9)
print "Making linear plot..."
for r in RValues:
    (zList_axis, fieldPICList_axis, fieldDirectList_axis) = extractRslice(r);
    plt.plot(zList_axis, fieldPICList_axis, 'r-', label="PIC field")
    plt.plot(zList_axis, fieldDirectList_axis, 'r--', label="Analytical field")
if useGrid:
    for iz in xrange(1,int(max(zList_axis))+1):
        print iz
        plt.axvline(iz,color='k', ls='--')
    plt.axvline(zWitness, color = 'g')
plt.xlabel("z [dz=%f um]" % (inputfile.dZ*1e4,))
plt.ylabel("Field [V/m]")

plt.figure(10)
print "Making log plot..."
for r in RValues:
    (zList_axis, fieldPICList_axis, fieldDirectList_axis) = extractRslice(r);
    plt.semilogy(zList_axis, fieldPICList_axis, 'r-')
    plt.semilogy(zList_axis, -fieldPICList_axis, 'b-')
    plt.semilogy(zList_axis, fieldDirectList_axis, 'r--')
    plt.semilogy(zList_axis, -fieldDirectList_axis, 'b--')
if useGrid:
    for iz in xrange(1,int(max(zList_axis))+1):
        plt.axvline(iz,color='k', ls='--')
    plt.axvline(zWitness, color = 'g')
plt.xlabel("z [dz=%f um]" % (inputfile.dZ*1e4,))
plt.ylabel("Field [V/m]")
print "done."

figureOffset = 11
order_last = fieldDirect_orderList[:,-1]
order_last = np.reshape(order_last,(-1, nrMax*pointsPerGrid))

order_fieldMax    = np.nanmax(fieldDirect_orderList);
order_fieldNegMax = np.nanmax(-fieldDirect_orderList);

fieldLevels_posMaxExp = int(np.log10(order_fieldMax))+1
fieldLevels_posSteps = fieldLevels_maxExp + 4 #starts at -3, all ints
fieldLevels_posSteps += fieldLevels_posSteps - 1  #add half-steps
fieldLevels_pos = np.logspace(-3, fieldLevels_posMaxExp, fieldLevels_posSteps);
print "FieldLevels_pos = "
print fieldLevels_pos
print

fieldLevels_negMaxExp = int(np.log10(order_fieldNegMax))+1
fieldLevels_negSteps = fieldLevels_negMaxExp + 4 #starts at -3, all ints incl. 0
fieldLevels_negSteps += fieldLevels_negSteps - 1    #add half-steps
fieldLevels_neg = np.logspace(-3, fieldLevels_negMaxExp, fieldLevels_negSteps);
print "FieldLevels_neg = "
print fieldLevels_neg
print

if mirrorOrders:
    for o in xrange(mirror_order+1):
        order_this = fieldDirect_orderList[:,o]
        order_this = np.reshape(order_this,(-1, nrMax*pointsPerGrid))
        order_relDiff = (order_this-order_last) / abs(order_last)
        
        plt.figure(figureOffset + o)
        CS1 = plt.contourf(zMat,rMat, order_this,
                           norm=LogNorm(), levels=fieldLevels_pos)
        plt.colorbar(CS1)
        CS2 = plt.contourf(zMat,rMat, -order_this, hatches='x',
                           norm=LogNorm(), levels=fieldLevels_neg)
        plt.colorbar(CS2)
        if hlineValue:
            plt.axhline(hlineValue, color='k')
        plt.xlabel("z [dz=%f um]" % (inputfile.dZ*1e4,))
        plt.ylabel("r [dz=%f um]" % (inputfile.dZ*1e4,))
        plt.title("Analytical field, order " + str(o));
        
        if o == mirror_order:
            continue;
        
        plt.figure(figureOffset+o+mirror_order+2)
        plt.contourf(zMat,rMat, np.abs(order_relDiff), norm=LogNorm())
        plt.colorbar()
        if hlineValue:
            plt.axhline(hlineValue, color='k')
        plt.xlabel("z [dz=%f um]" % (inputfile.dZ*1e4,))
        plt.ylabel("r [dz=%f um]" % (inputfile.dZ*1e4,))
        plt.title("abs(order %i - order %i) / abs(order %i)" % (o, mirror_order, mirror_order))
        
plt.show()
