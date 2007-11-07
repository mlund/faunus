#!/bin/bash

#
# Help script for the prot_dh.C protein-protein
# simulation program. Users should NOT modify
# this script! Is meant as a help to create
# automated scripts.
#
# M. Lund, Sep 2005
#

function setjobid() {
  jobid=$1
  infile="prot.inp${jobid}"
}

function newjob() {
  setjobid `date +".%k%M%S"`
}

#
#set defaults...
#
macro=10
randomseed=13
cell_r=100
maxsep=100
minsep=0
temp=298
dielec=78.69
prot_dp=7
prot_rot=2
hydrophob=0
penalty=0.
rotate="yes"
translate="yes"
minsnapshot="no"
hamaker=0
exe="${biosim}/src/prot_dh"
debye=900

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
randomseed   $randomseed
jobid        $jobid 
protein1     $protein1 
protein2     $protein2 
prot_dp      $prot_dp
prot_rot     $prot_rot
maxsep       $maxsep
minsep       $minsep
hamaker      $hamaker 
hydrophob    $hydrophob
rotate       $rotate 
translate    $translate
minsnapshot  $minsnapshot
macro        $macro
micro        $micro
debye        $debye
penalty      $penalty
" > $infile
}

