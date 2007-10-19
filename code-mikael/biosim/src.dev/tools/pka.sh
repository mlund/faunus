#!/bin/bash

function setjobid() {
  jobid=$1
  infile="pka.inp${jobid}"
}

function newjob() {
  setjobid `date +".%k%M%S"`
}

#
#set defaults...
#
r_dielec=0
macro=10
randomseed=13
cell_r=100
temp=298
dielec=78.69
dielec_i=78.69
ion_dp=50
nion1=0
chion1=1
nion2=0
chion2=-1
nion3=0
chion3=0
springconst=0.14    # = lB / (2*r_eq^3)
springeqdist=5
monomer_dp=3.
titrate="yes"
adjust_dp="yes"
ph=7
exe="${biosim}/src/pka"
springconst=0.5
movesites="no"
movecharges="no"
adjust_dp="yes"
watchsite=-1

newjob

#------------------------------------------------

function seed() {
  randomseed=`date +"%s"`
}

function simulate() {
  makeinput
  $exe $infile
};

function cleanup() {
  rm -fR .coord${jobid} $infile
}

function setsalt() {
  debye=`python -c "from math import sqrt; print '%.1f' % (3.04/sqrt($1))"`
}

function makeinput() {
echo "
cell_r       $cell_r
temp         $temp
dielec       $dielec
dielec_i     $dielec_i
r_dielec     $r_dielec
pH           $ph
randomseed   $randomseed
jobid        $jobid 
nion1        $nion1
nion2        $nion2 
nion3        $nion3 
chion1       $chion1 
chion2       $chion2 
chion3       $chion3 
ion_dp       $ion_dp
protein1     $protein1 
springconst  $springconst 
springeqdist $springeqdist
monomer_dp   $monomer_dp 
macro        $macro
micro        $micro
adjust_dp    $adjust_dp
springconst  $springconst
movesites    $movesites
movecharges  $movecharges
adjust_dp    $adjust_dp
watchsite    $watchsite
" > $infile
}

