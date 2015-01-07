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

if [ "$#" -ne 0 ] ; then
    echo "Usage: 2Dpic_restart.sh" > /dev/stderr
    echo "Generates total current for restarted 2D Arc-PIC runs. " > /dev/stderr
    echo "Run in folder run_name/gle " > /dev/stderr
    echo "Eg. 2Dpic_restart.sh"
    exit
fi


### RESCALING DATA TO DIMENSIONAL VALUES ###

dz=`grep "Grid size dz in Debyes:" ../input.txt | awk '{print $6}'`
dt=`grep "Timestep dt in Omega_pe-s:" ../input.txt | awk '{print $5}'`
Cext=`grep "External capacitance in F:" ../input.txt | awk '{print $5}'`
UNz=`grep "Potential at Zmax, UNz in Te:" ../input.txt | awk '{print $7}'`
U0=`grep "Potential at Zmin, U0 in Te:" ../input.txt | awk '{print $7}'`
echo "Cext= $Cext"
echo "U0= $U0 UNz= $UNz"

grep "n_ref:" ../readme.out | awk '{print $3}' > tmp.nref
grep "T_ref:" ../readme.out | awk '{print $3}' > tmp.Tref
grep "Ndb:" ../readme.out | awk '{print $3}' > tmp.Nd
grep "nav_dt:" ../readme.out | awk '{print $3}' > tmp.dt
grep "nstepsmin:" ../readme.out | awk '{print $3}' > tmp.min
last=`tail -1 ../IonWall1.dat | awk '{print $1}'`
paste tmp.nref tmp.Tref tmp.Nd tmp.dt tmp.min | awk '{printf "%.8e %.8e %.8e %d %d \n",$1,$2,$3,$4,$5;}' > restart.dat
rm tmp.nref tmp.Tref tmp.Nd tmp.dt tmp.min
runs=`wc -l restart.dat | awk '{print $1}'`




rm UI.dat out_counts.dat inw_counts.dat FE.dat
touch UI.dat out_counts.dat inw_counts.dat FE.dat

for ((i=1; i<$runs; i=i+1)); do

k=`echo $(($i+1))`

ne=`head -$i restart.dat | tail -1 | awk '{print $1}'`
Te=`head -$i restart.dat | tail -1 | awk '{print $2}'`
Nd=`head -$i restart.dat | tail -1 | awk '{print $3}'`
step=`head -$i restart.dat | tail -1 | awk '{print $4}'`
min=`head -$i restart.dat | tail -1 | awk '{print $5}'`
max=`head -$k restart.dat | tail -1 | awk '{printf "%d",$5-1;}'`

echo "Parameters run $i"
echo "ne= $ne Te= $Te"
echo "min= $min max= $max"
echo "step= $step Nd= $Nd"
echo "dz= $dz dt=$dt"

# OUTWARD PARTICLE COUNT

paste ../IonWall1.dat ../IonWall2.dat ../IonWall3.dat | awk -v max=$max -v min=$min -v dt=$dt -v ne=$ne '{omega=56414.6*sqrt(ne); t=$1; e1=$2; Cup1=$3; Cu1=$4; e2=$11; Cup2=$12; Cu2=$13; e3=$20; Cup3=$21; Cu3=$22; time=t*dt*1e+9/omega; if (($1<=max) && ($1>min)) {printf "%.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e \n", time,e1,Cup1,Cu1,e2,Cup2,Cu2,e3,Cup3,Cu3;}}' >> out_counts.dat


# INWARD PARTICLE COUNT

cat ../IonWall3.dat | awk -v max=$max -v min=$min -v dt=$dt -v ne=$ne '{omega=56414.6*sqrt(ne); t=$1; eFE=$5; eSEY=$6; CuY1=$7; CuY2=$8; Cuev=$9; time=t*dt*1e+9/omega; e=e*j*1e-8; Cu=Cu*j*1e-8; if (($1<=max) && ($1>min)) {printf "%.4e %.4e %.4e %.4e %.4e %.4e \n", time,eFE,eSEY,CuY1,CuY2,Cuev;}}' >> inw_counts.dat


# BETA EROSION, LOCAL FIELD (in GV/m)

cat ../IonWall1.dat | awk -v max=$max -v min=$min -v dt=$dt -v ne=$ne '{omega=56414.6*sqrt(ne); t=$1; B=$7; F=$8; time=t*dt*1e+9/omega; if (($1<=max) && ($1>min)) {printf "%.4e %.4e %.4e \n", time,B,F;}}' >> FE.dat


# TOTAL CURRENT (in A), EXTERNAL POTENTIAL (in kV)
# check current either through anode or cathode (better)

paste ../IonWall1.dat ../IonWall2.dat | awk -v max=$max -v min=$min -v D=$step -v dt=$dt -v ne=$ne -v Te=$Te -v Nd=$Nd '{dt_out=D*dt; I=2.693026e5*dt_out*Nd/Te/sqrt(Te); Qfac=1.519260e10*Nd*sqrt(ne)/Te/sqrt(Te); omega=56414.6*sqrt(ne); t=$1; e1=$2; Cup1=$3; ein1=$5; e2=$11; Cup2=$12; u=$16; Q=$17; time=t*dt*1e+9/omega; tot1=Cup1-e1+ein1; tot2=e2-Cup2; u=u*Te/1000; if (($1<=max) && ($1>min)) {printf "%.6e %.6e %.6e %.6e %.6e\n", time,tot1/I,tot2/I,u,Q/Qfac;}}' >> UI.dat


done

ne=`tail -1 restart.dat | awk '{print $1}'`
Te=`tail -1 restart.dat | awk '{print $2}'`
Nd=`tail -1 restart.dat | awk '{print $3}'`
step=`tail -1 restart.dat | awk '{print $4}'`
min=`tail -1 restart.dat | awk '{print $5}'`

echo "Parameters run $i"
echo "ne= $ne Te= $Te"
echo "min= $min max= $last"
echo "step= $step Nd= $Nd"
echo "dz= $dz dt=$dt"


# OUTWARD PARTICLE COUNT

paste ../IonWall1.dat ../IonWall2.dat ../IonWall3.dat | awk -v max=$last -v min=$min -v dt=$dt -v ne=$ne '{omega=56414.6*sqrt(ne); t=$1; e1=$2; Cup1=$3; Cu1=$4; e2=$11; Cup2=$12; Cu2=$13; e3=$20; Cup3=$21; Cu3=$22; time=t*dt*1e+9/omega; if (($1<=max) && ($1>min)) {printf "%.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e \n", time,e1,Cup1,Cu1,e2,Cup2,Cu2,e3,Cup3,Cu3;}}' >> out_counts.dat


# INWARD PARTICLE COUNT

cat ../IonWall3.dat | awk -v max=$last -v min=$min -v dt=$dt -v ne=$ne '{omega=56414.6*sqrt(ne); t=$1; eFE=$5; eSEY=$6; CuY1=$7; CuY2=$8; Cuev=$9; time=t*dt*1e+9/omega; e=e*j*1e-8; Cu=Cu*j*1e-8; if (($1<=max) && ($1>min)) {printf "%.4e %.4e %.4e %.4e %.4e %.4e \n", time,eFE,eSEY,CuY1,CuY2,Cuev;}}' >> inw_counts.dat


# BETA EROSION, LOCAL FIELD (in GV/m)

cat ../IonWall1.dat | awk -v max=$last -v min=$min -v dt=$dt -v ne=$ne '{omega=56414.6*sqrt(ne); t=$1; B=$7; F=$8;time=t*dt*1e+9/omega; if (($1<=max) && ($1>min)) {printf "%.4e %.4e %.4e \n", time,B,F;}}' >> FE.dat


# TOTAL CURRENT (in A), EXTERNAL POTENTIAL (in kV)
# check current either through anode or cathode (better)

paste ../IonWall1.dat ../IonWall2.dat | awk -v max=$last -v min=$min -v D=$step -v dt=$dt -v ne=$ne -v Te=$Te -v Nd=$Nd '{dt_out=D*dt; I=2.693026e5*dt_out*Nd/Te/sqrt(Te); Qfac=1.519260e10*Nd*sqrt(ne)/Te/sqrt(Te); omega=56414.6*sqrt(ne); t=$1; e1=$2; Cup1=$3; ein1=$5; e2=$11; Cup2=$12; u=$16; Q=$17; time=t*dt*1e+9/omega; tot1=Cup1-e1+ein1; tot2=e2-Cup2; u=u*Te/1000; if (($1<=max) && ($1>min)) {printf "%.6e %.6e %.6e %.6e %.6e\n", time,tot1/I,tot2/I,u,Q/Qfac;}}' >> UI.dat



# POWER CONSUMED

cat UI.dat | awk '{I=$2; u=$4; printf "%.4e \n", I*u*1000;}' > power_sign.dat
cat UI.dat | awk -v C=$Cext '{t=$1; I=$2; Q=$5; printf "%.4e %.4e \n",t,I*Q/C;}' > charge.dat

  #abs value for power
  awk '{ for (i=1; i<=NF; i=i+1) {if ($i<0) {$i=-$i;} print; }}' < power_sign.dat > power_sign2.dat
  paste UI.dat power_sign2.dat | awk '{t=$1; p=$6; printf "%.4e %.4e \n", t,p;}' > power_temp.dat

  lines=`wc -l power_temp.dat | awk '{print $1}'`
  lines=`echo $(($lines-1))`

  head -$lines power_temp.dat > power_1.dat
  tail -$lines power_temp.dat > power_2.dat
  paste power_1.dat power_2.dat | awk '{ t1=$1; p1=$2; t2=$3; p2=$4; P=(p2+p1)*(t2-t1)/2; E+=P; printf "%.4e %.4e %.4e \n",t2,P,1e-9*E;}' > power_integral.dat
  tail -1 power_integral.dat | awk '{printf "%s %.3e %s \n", "!energy",$3,"J";}' > power_header.dat
  cat power_header.dat power_temp.dat > power.dat

  head -$lines charge_temp.dat > charge_1.dat
  tail -$lines charge_temp.dat > charge_2.dat
  paste charge_1.dat charge_2.dat | awk '{ t1=$1; p1=$2; t2=$3; p2=$4; P=(p2+p1)*(t2-t1)/2; E+=P; printf "%.4e %.4e %.4e \n",t2,P,1e-9*E;}' > charge_integral.dat
  tail -1 charge_integral.dat | awk '{printf "%s %.3e %s \n", "!energy",$3,"J";}' > charge_header.dat
  cat charge_header.dat charge_temp.dat > charge.dat

  tail -1 power_header.dat
  tail -1 charge_header.dat
  writeE=`tail -1 power_header.dat | awk '{print $2}'`
  writeE2=`tail -1 charge_header.dat | awk '{print $2}'`

  eval "sed -e '/!write_here/ c\write \"$writeE J\" !write_here ' -e '/!write_2here/ c\write \"$writeE2 J\" !write_2here ' <power.gle >power_tmp.gle "
  rm power_1.dat power_2.dat power_header.dat power_temp.dat power_sign.dat power_sign2.dat
  rm charge_1.dat charge_2.dat charge_header.dat charge_temp.dat


gle -d png -o plot_UI.png UI.gle
gle -d png -o plot_FE.png FE.gle
gle -d png -o plot_P.png power_tmp.gle
gle -d png -o plot_e.png count_e.gle
gle -d png -o plot_Cup.png count_Cup.gle
gle -d png -o plot_Cu.png count_Cu.gle

echo "Done!"
