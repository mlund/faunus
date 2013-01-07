#!/bin/bash

# Submit with sbatch command

#SBATCH -t 168:00:00
#SBATCH -J rennes

# - Get exclusive access to whole node (16 cores on alarik)
#---SBATCH --exclusive

# - Number of cores
#SBATCH -n 8

# - Number of nodes
#BATCH -N 1

if [ "$SNIC_RESOURCE" == "alarik" ]
then
  module add openmpi/1.6.2/gcc/4.6.2
  cd $SLURM_SUBMIT_DIR
fi

corenum=8

exe=../mlund-rennes

function mkinput() {
echo "
Process  HASP  ASP     4.0     $pH
Process  HCTR  CTR     2.6     $pH
Process  HGLU  GLU     4.4     $pH
Process  HHIS  HIS     6.3     $pH
Process  HNTR  NTR     7.5     $pH
Process  HTYR  TYR     9.6     $pH
Process  HLYS  LYS     10.4    $pH
Process  HCYS  CYS     10.8    $pH
Process  HARG  ARG     12.0    $pH
" > eqtit

# ----------------------------------
#   GENERATE ATOM PARAMETERS
# ----------------------------------
echo "
Atom  Na      +1     2.0    0.1    1       no
Atom  Cl      -1     2.0    0.1    1       no
Atom  I       -1     2.0    0.1    1       no
Atom  SCN     -1     2.0    0.1    1       no

Atom  ASP     -1     3.6    0.1    110     no
Atom  HASP     0     3.6    0.1    110     no
Atom  CTR     -1     2.0    0.1    16      no
Atom  HCTR     0     2.0    0.1    16      no
Atom  GLU     -1     3.8    0.1    122     no
Atom  HGLU     0     3.8    0.1    122     no
Atom  HIS      0     3.9    0.1    130     no
Atom  HHIS     1     3.9    0.1    130     no
Atom  NTR      0     2.0    0.1    14      no
Atom  HNTR     1     2.0    0.1    14      no
Atom  TYR     -1     4.1    0.1    154     no
Atom  HTYR     0     4.1    0.1    154     no
Atom  LYS      0     3.7    0.1    116     no
Atom  HLYS     1     3.7    0.1    116     no
Atom  CYS     -1     3.6    0.1    103     no 
Atom  HCYS     0     3.6    0.1    103     no 
Atom  ARG      0     4.0    0.1    144     no
Atom  HARG     1     4.0    0.1    144     no

Atom  ALA      0     3.1    0.1    66      yes
Atom  ILE      0     3.6    0.1    102     yes
Atom  LEU      0     3.6    0.1    102     yes
Atom  MET      0     3.8    0.1    122     yes
Atom  PHE      0     3.9    0.1    138     yes
Atom  PRO      0     3.4    0.1    90      yes
Atom  TRP      0     4.3    0.1    176     yes
Atom  VAL      0     3.4    0.1    90      yes
Atom  SER      0     3.3    0.1    82      no
Atom  THR      0     3.5    0.1    94      no
Atom  ASN      0     3.6    0.1    108     no
Atom  GLN      0     3.8    0.1    120     no
Atom  GLY      0     2.9    0.1    54      no
" > cluster.atoms

for proc in {0..7}
do
if [ "$proc" == "0" ]; then salt=0.005; fi 
if [ "$proc" == "1" ]; then salt=0.010; fi 
if [ "$proc" == "2" ]; then salt=0.020; fi 
if [ "$proc" == "3" ]; then salt=0.040; fi 
if [ "$proc" == "4" ]; then salt=0.060; fi 
if [ "$proc" == "5" ]; then salt=0.080; fi 
if [ "$proc" == "6" ]; then salt=0.100; fi 
if [ "$proc" == "7" ]; then salt=0.120; fi 
echo "
atomlist               cluster.atoms
eq_processfile         eqtit
loop_macrosteps        10
loop_microsteps        $micro
cylinder_radius        $cylinder_radius
cylinder_len           $cylinder_len
temperature            298     # Kelvin
epsilon_r              78.3    # Water dielectric const
lj_eps                 0.05    # kT
polymer1_N             1
polymer2_N             1
polymer1_file          ../lfapo.aam
polymer2_file          ../blac.aam
transrot_transdp       $transdp
transrot_rotdp         2
dh_ionicstrength       $salt

temper_runfraction      $temper     # Set to one/zero to turn on/off tempering
temper_format           XYZQ 
" > mpi$proc.input
done
}

#rm -f mpi*

cylinder_len=350
cylinder_radius=80
pH=7
micro=1000
transdp=100
temper=0
mkinput
#mpiexec -np 8 $exe

temper=0.05
micro=2000000
transdp=20 #30
mkinput
mpiexec -np 8 $exe

