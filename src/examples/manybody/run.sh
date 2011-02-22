#!/bin/sh
#PBS -r n
#PBS -e q.err1
#PBS -o q.log1
### Mail to user
#PBS -m ae
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=1
# This job's working directory
cd $PBS_O_WORKDIR
echo Workin directory is `pwd`

echo Execution starts at   `date`

#--- Write config file ---
function mkinput() {
echo "
macrosteps $macrosteps
microsteps $microsteps
boxlen     $boxlen
bjerrum    $bjerrum
LJeps      $LJeps
nion1      $nion1
nion2      $nion2
nion3      $nion3
tion1      $tion1
tion2      $tion2
tion3      $tion3
atomfile   $atomfile
springconstant  $springconstant
springeqdist    $springeqdist
dp_monomer      $dp_monomer
dp_salt         $dp_salt
pH              $pH
pKa_core        $pKa_core
MUs NA -11.1 CL -11.1
protein1 ${protein1}
protein2 ${protein2}
nprot1   ${nprot1}
nprot2   ${nprot2}
" > manybody.conf
}
. use_modules
module add intel

#--- Input parameters ---
macrosteps=10
microsteps=1000
boxlen=145
bjerrum=7.12
LJeps=0.2
pH=7.5
nion3=100
nion2=300
nion1=0
tion2="NA"
tion3="CL"
tion1="GHOST"
protein1="lysozyme.aam"
protein2="apoalphaLa.aam"
nprot1=20
nprot2=20
atomfile="faunatoms.dat"
springconstant=1.0
springeqdist=4.0
dp_monomer=2
dp_salt=50
pKa_core=10.5
#
salt=139
for salt in 201 #100 300 900 
do
for pH in 7.5
do
  rm gcgroup.conf confout.aam
  suffix="test10mM"
  microsteps=10     
  nion2=$[100+760]
  nion3=$[100]
  nion1=$[0]
  mkinput
  ./manybody > ${suffix}.eq
  microsteps=100 
  mkinput
  ./manybody > ${suffix}.out
  mv smeared.aam ${suffix}.aam
done
done
