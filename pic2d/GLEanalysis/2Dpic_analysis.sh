#! /bin/bash
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
# Developers: Helga Timko, Kyrre Sjobak, Lotta Mether
#

if [ "$#" -ne 3 ] ; then
    echo "Usage: 2Dpic_analysis.sh <mintime> <maxtime> <every nth frame to analyse>" > /dev/stderr
    echo "Generates a full analysis of a 2D Arc-PIC c run. " > /dev/stderr
    echo "Run in folder run_name/gle " > /dev/stderr
    echo "Eg. to analyse all frames: 2Dpic_analysis.sh 10 200010 1"
    exit
fi

min=$1
max=$2
nth=$3

./2Dpic_joint.sh $min $max 2 5 $nth
mv pngs/joint/ pngs/joint_inner/
echo "Analysis joint_inner done!"

./2Dpic_joint.sh $min $max 50 100 $nth
mv pngs/joint/ pngs/joint_outer/
echo "Analysis joint_outer done!"

./2Dpic_joint.sh $min $max 20 40 $nth
echo "Analysis joint done!"

./2Dpic_temper.sh $min $max 10 50 $nth
echo "Analysis temper done!"

./2Dpic_coord.sh $min $max $nth
echo "Analysis coord done!"

./2Dpic_zoomcoord.sh $min $max $nth
echo "Analysis zoomcoord done!"

./2Dpic_efield.sh $min $max 0.5 $nth
echo "Analysis efield done!"

./2Dpic_pot.sh $min $max 0.5 $nth
echo "Analysis pot done!"

./2Dpic_current.sh $min $max
echo "Analysis current done!"
