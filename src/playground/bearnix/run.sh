#!/bin/bash

function mkinput() {

# Write titration and atom property files:
echo '
{
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
    "H-Lys" : { "bound":"HLYS" , "free":"LYS" , "pKd":10.4 , "pX":'$pH' }
  },
 
  "atomlist" :
  {
    "Na"   :  { "eps":0.124, "q": 1, "r":1.9, "mw":22.99 }, // epsilon in kJ/mol!
    "Cl"   :  { "eps":0.124, "q":-1, "r":1.7, "mw":35.45 },
    "X"    :  { "eps":0.010, "q":-1, "r":1.7, "mw":1.0 },
    "ASP"  :  { "eps":0.124, "q":-1, "r":3.6, "mw":110 },
    "HASP" :  { "eps":0.124, "q":0,  "r":3.6, "mw":110 },
    "LASP" :  { "eps":0.124, "q":2,  "r":3.6, "mw":110 },
    "CTR"  :  { "eps":0.124, "q":-1, "r":2.0, "mw":16 },
    "HCTR" :  { "eps":0.124, "q":0,  "r":2.0, "mw":16 },
    "GLU"  :  { "eps":0.124, "q":-1, "r":3.8, "mw":122 },
    "HGLU" :  { "eps":0.124, "q":0,  "r":3.8, "mw":122 },
    "LGLU" :  { "eps":0.124, "q":2,  "r":3.8, "mw":122 },
    "HIS"  :  { "eps":0.124, "q":0,  "r":3.9, "mw":130 },
    "HHIS" :  { "eps":0.124, "q":1,  "r":3.9, "mw":130 },
    "NTR"  :  { "eps":0.124, "q":0,  "r":2.0, "mw":14 },
    "HNTR" :  { "eps":0.124, "q":1,  "r":2.0, "mw":14 },
    "TYR"  :  { "eps":0.124, "q":-1, "r":4.1, "mw":154 },
    "HTYR" :  { "eps":0.124, "q":0,  "r":4.1, "mw":154 },
    "LYS"  :  { "eps":0.124, "q":0,  "r":3.7, "mw":116 },
    "HLYS" :  { "eps":0.124, "q":1,  "r":3.7, "mw":116 },
    "CYS"  :  { "eps":0.124, "q":-1, "r":3.6, "mw":103 },
    "HCYS" :  { "eps":0.124, "q":0,  "r":3.6, "mw":103 },
    "ARG"  :  { "eps":0.124, "q":0,  "r":4.0, "mw":144 },
    "HARG" :  { "eps":0.124, "q":1,  "r":4.0, "mw":144 },
    "ALA"  :  { "eps":0.124, "q":0,  "r":3.1, "mw":66,  "hydrophobic":true },
    "ILE"  :  { "eps":0.124, "q":0,  "r":3.6, "mw":102, "hydrophobic":true },
    "LEU"  :  { "eps":0.124, "q":0,  "r":3.6, "mw":102, "hydrophobic":true },
    "MET"  :  { "eps":0.124, "q":0,  "r":3.8, "mw":122, "hydrophobic":true },
    "PHE"  :  { "eps":0.124, "q":0,  "r":3.9, "mw":138, "hydrophobic":true },
    "PRO"  :  { "eps":0.124, "q":0,  "r":3.4, "mw":90,  "hydrophobic":true },
    "TRP"  :  { "eps":0.124, "q":0,  "r":4.3, "mw":176, "hydrophobic":true },
    "VAL"  :  { "eps":0.124, "q":0,  "r":3.4, "mw":90,  "hydrophobic":true },
    "SER"  :  { "eps":0.124, "q":0,  "r":3.3, "mw":82 },
    "THR"  :  { "eps":0.124, "q":0,  "r":3.5, "mw":94 },
    "ASN"  :  { "eps":0.124, "q":0,  "r":3.6, "mw":108 },
    "GLN"  :  { "eps":0.124, "q":0,  "r":3.8, "mw":120 },
    "GLY"  :  { "eps":0.124, "q":0,  "r":2.9, "mw":54 }
  }
}
' > cyl.json

# Write main input file:
echo "
atomlist               cyl.json
eq_processfile         cyl.json

loop_macrosteps        10
loop_microsteps        $micro

temperature            298     # Kelvin
epsilon_r              78.7    # Water dielectric const
dh_ionicstrength       $salt   # mol/l
eps_hydrophobic        0.5     # hydrophobic-hydrophobic LJ (kT)

cylinder_radius        $cylinder_radius # angstrom
cylinder_len           $cylinder_len    # angstrom

npt_P                  0       # mM
npt_dV                 0       # log(dV)
npt_runfraction        0.0
transrot_transdp       100     # Molecular translation parameter
transrot_rotdp         3       # Molecular rotation parameter
swapmv_runfraction     0.1     # Chance of performing titration

# Molecular species - currently only two different kinds
molecule1_N            1
molecule1              monopole.aam
molecule2_N            1
molecule2              dipole.aam
molecule_plane         1

movie                  yes     # save trajectory? (gromacs xtc format)
multipoledist_dr       0.2     # resolution of multipole analysis (angstrom)

# Atomic species - add up to ten different kinds
tion1                  Na
nion1                  0
dpion1                 10

tion2                  Cl
nion2                  0
dpion2                 10

" > cyl.input

# Generate some simple molecules:
# format: name anything x y z charge weight radius)
echo "2
 HIS  0   0.00   0.00   0.00    -0   1  3.0
 GLY  0   2.00   0.00   1.00     0   1  2.0
" > titratable.aam

echo "1
 X  0   0.00   0.00   0.00    +1   1  2.0
" > monopole.aam

echo "2
 X  0  -0.50   0.00   0.00    -1   1  2.0
 X  0   0.50   0.00   0.00     1   1  2.0
" > dipole.aam

echo "4
 X  0   1.00   1.00   0.00    +1   1  2.0
 X  0  -1.00  -1.00   0.00    +1   1  2.0
 X  0   1.00  -1.00   0.00    -1   1  2.0
 X  0  -1.00   1.00   0.00    -1   1  2.0
" > quadrupole.aam

}

exe=./bearnix-cyl
cylinder_len=50
cylinder_radius=60

# control equilibration and production runs
eqrun=true
prodrun=false
copy=true

for salt in 0.10  # loop over salt
do
  for pH in 7.0   # loop over pH
  do
    prefix="I${salt}-pH${pH}"

    if [ "$eqrun" = true ]; then
      echo "Equilibration run...(state file deleted)"
      rm -fR state
      micro=100000
      mkinput
      $exe
    fi

    if [ "$prodrun" = true ]; then
      echo "Production run..."
      micro=100000
      mkinput
      $exe > $prefix.out
    fi

    if [ "$copy" = true ]; then
      cp "$0" $prefix.sh
      cp state $prefix.state
      cp traj.xtc $prefix.xtc
      cp cyl.input $prefix.input
      cp rdf_p2p.dat $prefix.rdf 
      cp confout.pqr $prefix.pqr
      cp multipole.dat $prefix.multipole
    fi
  done
done
