#!/bin/bash

exit "not finished..."
exit

# .atom -> .json
# gawk '{print "\""$2"\" :  { \"q\":"$3",  \"r\":"$4", \"mw\":"$6" },"}'

#
# Interactive script to titrate a rigid molecular structures
# using an implicit salt description.
#

faunus=$HOME/faunus/
pkaexe=src/tools/titrate_implicit
base=pka

function mktit() {
echo "
Process  H3PO4 H2PO4   2.12    $pH
Process  H2PO4 HPO4    7.21    $pH
Process  HPO4  PO4     12.67   $pH
Process  HNTR  NTR     7.5     $pH
Process  HTYR  TYR     9.6     $pH
Process  HLYS  LYS     10.4    $pH
Process  HCYS  CYS     10.8    $pH
" > $base.titration
}

# ----------------------------------
#   GENERATE ATOM PARAMETERS
# ----------------------------------
function mkatoms() {
echo "
{
  "processes" :
  {
    "H-Asp" : { "bound":"HASP" , "free":"ASP" , "pKd":4.0  , "pKx":$pH }
    "H-Ctr" : { "bound":"HCTR" , "free":"CTR" , "pKd":2.6  , "pKx":$pH }
    "H-Glu" : { "bound":"HGLU" , "free":"GLU" , "pKd":4.4  , "pKx":$pH }
    "H-His" : { "bound":"HHIS" , "free":"HIS" , "pKd":6.3  , "pKx":$pH }
    "H-Arg" : { "bound":"HARG" , "free":"ARG" , "pKd":12.0 , "pKx":$pH }
  }

  "atomlist" :
  {
    "Na"   :  { "q": 1, "r":1.9, "mw":22.99 }, // sodium ion
    "Cl"   :  { "q":-1, "r":1.7, "mw":35.45 }, // chloride ion
    "I"    :  { "q":-1, "r":2.0, "mw":1 },
    "SCN"  :  { "q":-1, "r":2.0, "mw":1 },
    "ASP"  :  { "q":-1, "r":3.6, "mw":110 },
    "HASP" :  { "q":0,  "r":3.6, "mw":110 },
    "CTR"  :  { "q":-1, "r":2.0, "mw":16 },
    "HCTR" :  { "q":0,  "r":2.0, "mw":16 },
    "GLU"  :  { "q":-1, "r":3.8, "mw":122 },
    "HGLU" :  { "q":0,  "r":3.8, "mw":122 },
    "HIS"  :  { "q":0,  "r":3.9, "mw":130 },
    "HHIS" :  { "q":1,  "r":3.9, "mw":130 },
    "NTR"  :  { "q":0,  "r":2.0, "mw":14 },
    "HNTR" :  { "q":1,  "r":2.0, "mw":14 },
    "TYR"  :  { "q":-1, "r":4.1, "mw":154 },
    "HTYR" :  { "q":0,  "r":4.1, "mw":154 },
    "LYS"  :  { "q":0,  "r":3.7, "mw":116 },
    "HLYS" :  { "q":1,  "r":3.7, "mw":116 },
    "CYS"  :  { "q":-1, "r":3.6, "mw":103 },
    "HCYS" :  { "q":0,  "r":3.6, "mw":103 },
    "ARG"  :  { "q":0,  "r":4.0, "mw":144 },
    "HARG" :  { "q":1,  "r":4.0, "mw":144 },
    "ALA"  :  { "q":0,  "r":3.1, "mw":66 },
    "ILE"  :  { "q":0,  "r":3.6, "mw":102 },
    "LEU"  :  { "q":0,  "r":3.6, "mw":102 },
    "MET"  :  { "q":0,  "r":3.8, "mw":122 },
    "PHE"  :  { "q":0,  "r":3.9, "mw":138 },
    "PRO"  :  { "q":0,  "r":3.4, "mw":90 },
    "TRP"  :  { "q":0,  "r":4.3, "mw":176 },
    "VAL"  :  { "q":0,  "r":3.4, "mw":90 },
    "SER"  :  { "q":0,  "r":3.3, "mw":82 },
    "THR"  :  { "q":0,  "r":3.5, "mw":94 },
    "ASN"  :  { "q":0,  "r":3.6, "mw":108 },
    "GLN"  :  { "q":0,  "r":3.8, "mw":120 },
    "GLY"  :  { "q":0,  "r":2.9, "mw":54 }
  }
}

{ "atomnames"
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
  " > $base.atoms
}

function mkinput() {
echo "
atomlist               $base.atoms
eq_processfile         $base.titration
loop_macrosteps        10
loop_microsteps        $micro
temperature            298     # Kelvin
epsilon_r              78.7    # Water dielectric const
lj_eps                 0.0     # kT
dh_ionicstrength       $ionicstr   # mol/l
sphere_radius          500     # Angstroms
molecule               $struct
" > input
}

function userinput() {
read -p "Path to faunus directory ["$faunus"]: " i
faunus=${i:-$faunus}
if [ ! -d $faunus/src/tools ]; then
  echo "Invalid faunus directory." ; exit 1
fi
if [ ! -f $faunus/$pkaexe ]; then
  echo "Could not find titration program. Did you remember to make it?" ; exit 1
fi

read -p "AAM structure file to titrate: " struct
if [ ! -f $struct ]; then
  echo "Structure file does not exist." ; exit 1
fi

ionicstr=0.020
read -p "Ionic strength ["$ionicstr" mol/l]: " i
ionicstr=${i:-$ionicstr}

pHrange="2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0"
read -p "pH range [$pHrange]: " i
pHrange=${i:-$pHrange}

echo
echo "Summary of simulation:"
echo "  Titration program        = $faunus/$pkaexe"
echo "  Molecular structure file = $struct"
echo "  Ionic strength           = $ionicstr mol/l"
echo "  pH range                 = $pHrange"
}

userinput

echo
read -p "Start simulations? "
echo

echo "pH    Z        C       mu"
for pH in $pHrange; do
  echo -n "$pH "
  mkatoms
  mktit

  # Equilibration run
  rm -f state
  micro=1000
  mkinput
  $faunus/$pkaexe > eq

  # Production run
  micro=10000
  prefix="pH$pH-salt$ionicstr"
  mkinput
  $faunus/$pkaexe > $prefix.out
  mv confout.pqr $prefix.pqr
  Z=`grep "Analysis: Multipole" -A 6 $prefix.out | tail -n 1  | gawk '{print $2}'`
  C=`grep "Analysis: Multipole" -A 6 $prefix.out | tail -n 1  | gawk '{print $3}'`
  mu=`grep "Analysis: Multipole" -A 6 $prefix.out | tail -n 1  | gawk '{print $4}'`
  echo $Z $C $mu
done

