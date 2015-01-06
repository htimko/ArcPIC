#!/usr/bin/env python
#
# Copyright 2010-2014 CERN and Helsinki Institute of Physics.
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

print "************************************************"
print "*** makeRundir.py:                           ***"
print "*** Small utility for copying the sources    ***"
print "*** and standard input file to a new path.   ***"
print "*** Usage: makeRundir.py new/path/           ***"
print "************************************************"

import sys, os, shutil, stat
import subprocess
import datetime

newpath = sys.argv[1]
if newpath[-1] == "/":
    newpath = newpath[:-1]
(newpathroot, runname) = os.path.split(newpath)

if os.path.exists(newpath):
    print "Path \""+newpath+"\" already exists"
    exit(1)
if not os.path.isdir(newpathroot):
    print "Path \""+newpathroot+"\" does not exist."
    exit(1)

print "Copying files..."
os.mkdir(newpath)
shutil.copy("input.txt", os.path.join(newpath,"input.txt"))
shutil.copy("calcScaling.py", os.path.join(newpath,"calcScaling.py"))
def ignoFunc(src,dirlist):
    rl = []
    for l in dirlist:
        if not (l.endswith(".cpp") or l.endswith(".h") or l.endswith(".py") or l.endswith(".dat")):
            rl.append(l);
    return rl
shutil.copytree("h", os.path.join(newpath,"h"), ignore=ignoFunc)
shutil.copytree("src", os.path.join(newpath,"src"), ignore=ignoFunc)
shutil.copytree("tests", os.path.join(newpath,"tests"), ignore=ignoFunc)
shutil.copytree("inputdata", os.path.join(newpath,"inputdata"), ignore=ignoFunc)

print "Modifying the Makefile with new $RUNNAME and $WORK..."
mkfile = open("src/Makefile",'r')
mkfile2 = open(os.path.join(newpath,"src","Makefile"),'w')
for line in mkfile:
    l = line
    if "MAKERUNDIR" in line:
        if line.startswith("RUNNAME"):
            l = "RUNNAME = " + runname + "\n"
        elif line.startswith("WORK"):
            l = "WORK = " + newpath + "/\n"
        elif line.startswith("SUPPORTroot"):
            arcpicroot = os.path.split(os.path.split(os.getcwd())[0])[0] #two dirs up
            l = "SUPPORTroot = " + os.path.join(arcpicroot, "support", "pic2d", "trunk") + "\n"
        else:
            print "Magic comment in the wrong place"
    mkfile2.write(l);
mkfile.close()
mkfile2.close()

print "Installing GLE analysis scripts..."
SVNbasedir = os.path.join("..","..") #Change this if using a branch! (should in principle be deductible from 'svn info' output)
anaScriptsFolder = os.path.join(SVNbasedir, "analysis", "pic2d", "trunk")
anaSupportFolder = os.path.join(SVNbasedir, "support", "analysis", "pic2d", "gle_template_redblue", "trunk")
os.mkdir(os.path.join(newpath, "GLEanalysis"))

def copyFolderContents(srcFolder, targetFolder):
    "Adapted from shutil.copytree"
    names = os.listdir(srcFolder)
    for name in names:
        srcname = os.path.join(srcFolder, name)
        targetname = os.path.join(targetFolder, name)
        if os.path.isdir(srcname):
            if name == ".svn":
                continue #Skip these
            #Makes directory and copys the stuff
            shutil.copytree(srcname,targetname)
        else:
            shutil.copy2(srcname,targetname)

# COMMENT OUT FOR NOW
# copyFolderContents(anaScriptsFolder, os.path.join(newpath, "GLEanalysis"))

#Chmod every file which has a hasbang to u+x
for name in os.listdir(os.path.join(newpath, "GLEanalysis")):
    fn = os.path.join(newpath, "GLEanalysis",name)
    if not os.path.isfile(fn):
        continue
    
    fno = open(fn,'r');
    fnoline = fno.readline();
    fno.close()
    if not fnoline.startswith("#!"):
        continue

    os.chmod(fn, os.stat(fn).st_mode | stat.S_IXUSR)

#copyFolderContents(anaSupportFolder, os.path.join(newpath, "GLEanalysis"))

#Make a SVN status log
#print "Making SVN status log..."
#svnlog = open(os.path.join(newpath, "SVNstatus.txt"), 'w');
#svnlog.write("SVN log created by makerundir at " + str(datetime.datetime.now()) + "\n")
#svnlog.write("Output from 'svn info':\n")
#svnlog.write("#########################################\n")
#svnlog.flush()
#subprocess.call(["svn", "info"], stdout=svnlog, stderr=subprocess.STDOUT)
#svnlog.write("Output from 'svnversion':\n")
#svnlog.write("#########################################\n")
#svnlog.flush()
#subprocess.call(["svnversion"], stdout=svnlog, stderr=subprocess.STDOUT)
#svnlog.write("#########################################\n")
#svnlog.write("Output from 'svn diff':\n")
#svnlog.write("#########################################\n")
#svnlog.flush()
#subprocess.call(["svn", "diff"], stdout=svnlog, stderr=subprocess.STDOUT)
#svnlog.write("#########################################\n")

print "Installed ArcPic in \"" + newpath + "\""
