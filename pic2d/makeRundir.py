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
        if not (l.endswith(".cpp") or l.endswith(".h") or l.endswith(".py") or \
                l.endswith(".dat") or l.endswith(".gle") or l.endswith(".sh")):
            rl.append(l);
    return rl
shutil.copytree("h", os.path.join(newpath,"h"), ignore=ignoFunc)
shutil.copytree("src", os.path.join(newpath,"src"), ignore=ignoFunc)
shutil.copytree("tests", os.path.join(newpath,"tests"), ignore=ignoFunc)
shutil.copytree("inputdata", os.path.join(newpath,"inputdata"), ignore=ignoFunc)
shutil.copytree("GLEanalysis", os.path.join(newpath,"GLEanalysis"), ignore=ignoFunc)
os.mkdir(os.path.join(newpath, "GLEanalysis", "pngs"))


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

#Make a GIT status log
print "Making GIT status log..."
gitlog = open(os.path.join(newpath, "GITstatus.txt"), 'w');
gitlog.write("GIT log created by makerundir at " + str(datetime.datetime.now()) + "\n")
gitlog.write("\n\n")

gitlog.write("Output from 'git log -n 1':\n")
gitlog.write("#########################################\n")
gitlog.flush()
subprocess.call(["git", "log", "-n", "1"], stdout=gitlog, stderr=subprocess.STDOUT)
gitlog.write("\n\n")

gitlog.write("Output from 'git status':\n")
gitlog.write("#########################################\n")
gitlog.flush()
subprocess.call(["git", "status"], stdout=gitlog, stderr=subprocess.STDOUT)
gitlog.write("\n\n")

gitlog.write("Output from 'git diff':\n")
gitlog.write("#########################################\n")
gitlog.flush()
subprocess.call(["git", "diff"], stdout=gitlog, stderr=subprocess.STDOUT)
gitlog.write("\n\n")

print "Installed ArcPic in \"" + newpath + "\""
