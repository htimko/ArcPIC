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

if [ "$#" -ne 5 ] ; then
    echo "Usage: 2Dpic_temper.sh <mintime> <maxtime> <cross-section NG2> <cross-section NG3> <every nth frame to analyse>" > /dev/stderr
    echo "Generates joint movies of T for 2D Arc-PIC code. " > /dev/stderr
    echo "Run in folder run_name/gle " > /dev/stderr
    echo "Eg. to analyse all frames: 2Dpic_temper.sh 10 200010 2 5 1"
    exit
fi

min=$1
max=$2
grid2=$3
grid3=$4
nth=$5

grid1=0

step1=`head -2 ../out/timeIndex.dat | tail -1 | awk '{print $1}'`
step2=`head -3 ../out/timeIndex.dat | tail -1 | awk '{print $1}'`
step=`echo | awk -v step1=$step1 -v step2=$step2 -v nth=$nth '{ printf "%d", step2-step1*nth}'`
code=`echo | awk -v m=$min '{ printf "%08ld", m}'`

echo "Starting code is $code"
echo "Timestep is $step"


### RESCALING DATA TO DIMENSIONAL VALUES ###

ne=`grep "Reference density in 1/cm3:" ../input.txt | awk '{print $5}'`
Te=`grep "Reference temperature in eV:" ../input.txt | awk '{print $5}'`
nr=`grep "Number of cells nr, nz:" ../input.txt | awk '{print $6}'`
nz=`grep "Number of cells nr, nz:" ../input.txt | awk '{print $7}'`
dz=`grep "Grid size dz in Debyes:" ../input.txt | awk '{print $6}'`
dt=`grep "Timestep dt in Omega_pe-s:" ../input.txt | awk '{print $5}'`
dt_out=`grep "Outputting time dt_out in O_pe-s:" ../input.txt | awk '{print $6}'`

echo "Parameters"
echo "ne= $ne Te= $Te"
echo "nr= $nr nz= $nz"
echo "dz= $dz dt= $dt"
echo "dt_out= $dt_out"

### *** ###


# RESCALE GRID TO DIMENSIONLESS DISTANCE
grd1=`echo | awk -v dz=$dz -v g=$grid1 '{ printf "%.5e", g*dz}'`
grd2=`echo | awk -v dz=$dz -v g=$grid2 '{ printf "%.5e", g*dz}'`
grd3=`echo | awk -v dz=$dz -v g=$grid3 '{ printf "%.5e", g*dz}'`

echo "grid1,2,3: $grid1 $grid2 $grid3"
echo "grd1,2,3: $grd1 $grd2 $grd3 "


# WRITING CORRECT DISTANCES INTO GLE
  echo | awk -v g1=$grd1 -v n=$ne -v Te=$Te '{Ld=sqrt(552635*Te/n); g1=g1*Ld*10000; printf"%2.1f \n",g1;}' > grid.tmp
  echo | awk -v g2=$grd2 -v n=$ne -v Te=$Te '{Ld=sqrt(552635*Te/n); g2=g2*Ld*10000; printf"%2.1f \n",g2;}' >> grid.tmp
  echo | awk -v g3=$grd3 -v n=$ne -v Te=$Te '{Ld=sqrt(552635*Te/n); g3=g3*Ld*10000; printf"%2.1f \n",g3;}' >> grid.tmp
    g1=`head -1 grid.tmp | awk '{print $1}'`
    g2=`head -2 grid.tmp | tail -1 | awk '{print $1}'`
    g3=`head -3 grid.tmp | tail -1 | awk '{print $1}'`
    rm grid.tmp
    echo "g1,2,3: $g1 $g2 $g3 "
    
    # echo "At sed grid"
    eval "sed -e '/grid1/ c\write \"$g1 um\" ' <temper.gle >temper_tmp.gle "
    eval "sed -e '/grid2/ c\write \"$g2 um\" ' <temper_tmp.gle >temper_tmp2.gle "
    eval "sed -e '/grid3/ c\write \"$g3 um\" ' <temper_tmp2.gle >temper_grid.gle "
    rm temper_tmp.gle temper_tmp2.gle


# PRINT OUT TIME (in ns)
rm time.dat
cat ../out/timeIndex.dat | awk -v dt=$dt -v ne=$ne '{omega=56414.6*sqrt(ne); t=$1; time=t*dt*1e+9/omega; printf "%08ld %1.2f \n", t, time}' > time.dat

rm -r pngs/temper
mkdir pngs/temper



k=1
for ((i=$min; i<=$max; i=i+$step)); do

  now_is=`grep "$code" time.dat | awk '{print $2}'`
  echo "now_is $now_is"
  echo "code is $code"
	
### TEMPERATURE Z-DIR ###

    rm Tez_a.dat Tez_b.dat Tez_c.dat Tiz_a.dat Tiz_b.dat Tiz_c.dat 

    # create temporary file with rescaled values, given grids in r-dir
    cat ../out/Tez${code}.dat | awk -v n=$ne -v Te=$Te -v g=$grd1 '{ Ld=sqrt(552635*Te/n); z=$2; T=$3; z=z*Ld*10000; T=T*Te; if($1==g) {printf "%.4f %.6e \n", z,T;}}' > Tez_a.dat
    cat ../out/Tez${code}.dat | awk -v n=$ne -v Te=$Te -v g=$grd2 '{ Ld=sqrt(552635*Te/n); z=$2; T=$3; z=z*Ld*10000; T=T*Te; if($1==g) {printf "%.4f %.6e \n", z,T;}}' > Tez_b.dat
    cat ../out/Tez${code}.dat | awk -v n=$ne -v Te=$Te -v g=$grd3 '{ Ld=sqrt(552635*Te/n); z=$2; T=$3; z=z*Ld*10000; T=T*Te; if($1==g) {printf "%.4f %.6e \n", z,T;}}' > Tez_c.dat
	
    cat ../out/Tiz${code}.dat | awk -v n=$ne -v Te=$Te -v g=$grd1 '{ Ld=sqrt(552635*Te/n); z=$2; T=$3; z=z*Ld*10000; T=T*Te; if($1==g) {printf "%.4f %.6e \n", z,T;}}' > Tiz_a.dat
    cat ../out/Tiz${code}.dat | awk -v n=$ne -v Te=$Te -v g=$grd2 '{ Ld=sqrt(552635*Te/n); z=$2; T=$3; z=z*Ld*10000; T=T*Te; if($1==g) {printf "%.4f %.6e \n", z,T;}}' > Tiz_b.dat
    cat ../out/Tiz${code}.dat | awk -v n=$ne -v Te=$Te -v g=$grd3 '{ Ld=sqrt(552635*Te/n); z=$2; T=$3; z=z*Ld*10000; T=T*Te; if($1==g) {printf "%.4f %.6e \n", z,T;}}' > Tiz_c.dat
	


### TOTAL TEMPERATURE ###
    
    rm Tetot_a.dat Tetot_b.dat Tetot_c.dat Titot_a.dat Titot_b.dat Titot_c.dat
    
    # create temporary file with rescaled values
    paste ../out/Tez${code}.dat ../out/Ter${code}.dat ../out/Tet${code}.dat | awk -v n=$ne -v Te=$Te -v g=$grd1 '{ Ld=sqrt(552635*Te/n); z=$2; T=$3+$6+$9; z=z*Ld*10000; T=T*Te; if($1==g) {printf "%.4f %.6e \n", z,T;}}' > Tetot_a.dat
    paste ../out/Tez${code}.dat ../out/Ter${code}.dat ../out/Tet${code}.dat | awk -v n=$ne -v Te=$Te -v g=$grd2 '{ Ld=sqrt(552635*Te/n); z=$2; T=$3+$6+$9; z=z*Ld*10000; T=T*Te; if($1==g) {printf "%.4f %.6e \n", z,T;}}' > Tetot_b.dat
    paste ../out/Tez${code}.dat ../out/Ter${code}.dat ../out/Tet${code}.dat | awk -v n=$ne -v Te=$Te -v g=$grd3 '{ Ld=sqrt(552635*Te/n); z=$2; T=$3+$6+$9; z=z*Ld*10000; T=T*Te; if($1==g) {printf "%.4f %.6e \n", z,T;}}' > Tetot_c.dat
	
    paste ../out/Tiz${code}.dat ../out/Tir${code}.dat ../out/Tit${code}.dat | awk -v n=$ne -v Te=$Te -v g=$grd1 '{ Ld=sqrt(552635*Te/n); z=$2; T=$3+$6+$9; z=z*Ld*10000; T=T*Te; if($1==g) {printf "%.4f %.6e \n", z,T;}}' > Titot_a.dat
    paste ../out/Tiz${code}.dat ../out/Tir${code}.dat ../out/Tit${code}.dat | awk -v n=$ne -v Te=$Te -v g=$grd2 '{ Ld=sqrt(552635*Te/n); z=$2; T=$3+$6+$9; z=z*Ld*10000; T=T*Te; if($1==g) {printf "%.4f %.6e \n", z,T;}}' > Titot_b.dat
    paste ../out/Tiz${code}.dat ../out/Tir${code}.dat ../out/Tit${code}.dat | awk -v n=$ne -v Te=$Te -v g=$grd3 '{ Ld=sqrt(552635*Te/n); z=$2; T=$3+$6+$9; z=z*Ld*10000; T=T*Te; if($1==g) {printf "%.4f %.6e \n", z,T;}}' > Titot_c.dat



		

### execute gle, into png or jpg format
    #echo "At sed time"
    eval "sed -e '/Time/ c\write \"Time $now_is ns\" ' <temper_grid.gle >temper_tmp.gle "


#   gle -d jpg -o joint${k} joint_tmp.gle
#   mv joint${k}.jpg pngs/joint/.

# HIGH RESOLUTION/LOW RESOLUTION
#   gle -d png -dpi 300 -o joint${k} joint_tmp.gle
    kk=`echo | awk -v k=$k '{ printf "%03ld", k}'`
    gle -d png -o temper${kk} temper_tmp.gle
    mv temper${kk}.png pngs/temper/.

    rm temper_tmp.gle

  k=`echo $(($k+1))`
  code=`echo | awk -v i=$i -v s=$step '{x=i+s; printf "%08ld", x;}'`
  echo " "
		


done

#rm joint_grid.gle
ffmpeg -sameq -r 20 -f image2 -i pngs/temper/temper%03d.png  pngs/temper/movie_temper.mpg

echo "Done with analysis!"
