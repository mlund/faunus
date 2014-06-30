#!/bin/bash

function mkatoms() {
echo '
{

  "partition coefficients" :
  {
    "Na" :
    {
      "ASP" : { "Kp" : 0.5 },
      "BB"  : { "Kp" : 1.5 }
    }
  }

  "processes" :
  {
    "H-Asp" : { "bound":"HASP" , "free":"ASP" , "pKd":4.0  , "pX":'$pH' },
    "H-Ctr" : { "bound":"HCTR" , "free":"CTR" , "pKd":2.6  , "pX":'$pH' },
    "H-Glu" : { "bound":"HGLU" , "free":"GLU" , "pKd":4.4  , "pX":'$pH' },
    "H-His" : { "bound":"HHIS" , "free":"HIS" , "pKd":6.3  , "pX":'$pH' },
    "H-Arg" : { "bound":"HARG" , "free":"ARG" , "pKd":12.0 , "pX":'$pH' },
    "H-Ntr" : { "bound":"HNTR" , "free":"NTR" , "pKd":7.5  , "pX":'$pH' },
    "H-Cys" : { "bound":"HCYS" , "free":"CYS" , "pKd":10.8 , "pX":'$pH' },
    "H-Tyr" : { "bound":"HTYR" , "free":"TYR" , "pKd":9.6  , "pX":'$pH' },
    "H-Lys" : { "bound":"HLYS" , "free":"LYS" , "pKd":10.4 , "pX":'$pH' },
    "K1"    : { "bound":"H3PO4", "free":"H2PO4","pKd":2.12,  "pX":'$pH' },
    "K2"    : { "bound":"H2PO4", "free":"HPO4", "pKd":7.21,  "pX":'$pH' },
    "K3"    : { "bound":"HPO4",  "free":"PO4",  "pKd":12.67, "pX":'$pH' }
  },
 
  "atomlist" :
  {
    "H3PO4":  { "q":0,  "r":2.0 },           // simple ions
    "H2PO4":  { "q":-1, "r":2.0 },
    "HPO4" :  { "q":-2, "r":2.0 },
    "PO4"  :  { "q":-3, "r":2.0 },
    "Na"   :  { "q": 1, "r":1.9, "mw":22.99 },
    "Cl"   :  { "q":-1, "r":1.7, "mw":35.45 },
    "I"    :  { "q":-1, "r":2.0, "mw":1 },
    "SCN"  :  { "q":-1, "r":2.0, "mw":1 },
    "ASP"  :  { "q":-1, "r":3.6, "mw":110 }, // amino acids
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
' > titrate.json

#----------------------------------------------------
# Make file with molecule of ONE atom
#----------------------------------------------------
echo "1
 PO4  0   0.00   0.00   0.00  -3.0   1.0  2.0
" > titrate.aam
}

function mkinput() {
echo "
atomlist               titrate.json
eq_processfile         titrate.json
molecule               titrate.aam
loop_macrosteps        10
loop_microsteps        10000
temperature            298     # Kelvin
epsilon_r              78.7    # Water dielectric const
lj_eps                 0.0     # kT
dh_ionicstrength       0.1     # mol/l
sphere_radius          500     # Angstroms
" > titrate.input
}

for pH in 7.21; do
  mkatoms
  rm -f state
  mkinput
  ./titrate
done

