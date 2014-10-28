#!/bin/bash

function mkinput() {

echo '
{
  "atomlist" :
  {
    "PROT" :  { "q": 0, "r":10., "mw":99999 },
    "MU"   :  { "q": 0, "r":0.0, "mw":0.001 },
    "Na"   :  { "q": 1, "r":2.0, "mw":1.000 },
    "Cl"   :  { "q":-1, "r":2.0, "mw":1.000 }
  }
}
' > cyl.json

# Write main input file:
echo "
atomlist               cyl.json

loop_macrosteps        10
loop_microsteps        $micro

temperature            298     # Kelvin
epsilon_r              78.7    # Water dielectric const

cylinder_radius        $cylinder_radius # angstrom
cylinder_len           $cylinder_len    # angstrom

transrot_transdp       10      # Molecular translation parameter
transrot_rotdp         3       # Molecular rotation parameter

# Molecules
molecule1_N            1
molecule1              mol.pqr
molecule2_N            1
molecule2              mol.pqr
muscalar               $muscalar # magnitude of off-center dipoles

molecule_plane         yes
movie                  no        # save trajectory? (gromacs xtc format)

tion1                  Na
nion1                  $Nion
tion2                  Cl
nion2                  $Nanion
mv_particle_genericdp  40

" > cyl.input

# Charged colloid with off-center dipole, 1 angstrom below surface
# atom# atomname resname res#    x        y        z      q     r
echo "
ATOM      1 PROT PROT    1       0.000    0.000    0.000  10.00 10.00
ATOM      2 MU   MU      2       0.000    0.000    9.000  0.000 0.000
" > mol.pqr
}

exe=./cyl
cylinder_len=200
cylinder_radius=30

# control equilibration and production runs
eqrun=true
prodrun=true
copy=true

for muscalar in 1.0 #2.0    # off-center dipole moment (Debye), 90 deg. hard coded
do
  for salt in 0.020 #0.040 0.080 0.160 # 1:1 salt concentraion (mol/l)
  do
    vol=`python -c "print 3.1416*($cylinder_radius)**2 * $cylinder_len"`
    conc=`python -c "print 1e-27 * 6.022e23 * $salt"`
    Nion=`python -c "print int(round($conc * $vol))"`
    prefix="salt${salt}-mu$muscalar"

    Nanion=`python -c "print $Nion+2*10"`  # number of chloride ions

    if [ "$eqrun" = true ]; then
      echo "Equilibration run...(state file deleted)"
      rm -fR state
      micro=2000
      mkinput
      $exe > $prefix.eq
    fi

    if [ "$prodrun" = true ]; then
      echo "Production run..."
      micro=1000000
      mkinput
      $exe > $prefix.out
    fi

    if [ "$copy" = true ]; then
      cp "$0" $prefix.sh
      cp state $prefix.state
      cp rdf_p2p.dat $prefix.rdf
      cp confout.pqr $prefix.pqr
    fi
  done
done
