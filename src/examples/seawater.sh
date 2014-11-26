#!/bin/bash

# THIS RUN SCRIPT IS USED AS A UNIT TEST SO PLEASE
# DO NOT UPLOAD ANY MODIFIED VERSIONS TO GIT UNLESS
# TO UPDATE THE TEST.

#source_tests_dir="`dirname $0`"
#cp -f $source_tests_dir/bulk.test . 2> /dev/null
#cp -f $source_tests_dir/bulk.state state 2> /dev/null

echo '{
  "atomlist" : {
    "Na"   : { "q": 1.0, "r":1.40, "dp":40. },
    "Cl"   : { "q":-1.0, "r":1.81, "dp":40. },
    "Mg"   : { "q": 2.0, "r":2.95, "dp":10. },
    "Ca"   : { "q": 2.0, "r":2.74, "dp":10. },
    "SO4"  : { "q":-2.0, "r":1.30, "dp":10. },
    "HCO3" : { "q":-1.0, "r":0.90, "dp":40. },
    "CO3"  : { "q":-2.0, "r":1.30, "dp":10. }
  }
}' > seawater.json

echo "
atomlist        seawater.json    # atom properties
cuboid_len      100          # angstrom

temperature     298.15       # K
epsilon_r       80           # dielectric const.
coulomb_cut     14.          # coulomb cutoff [angstrom]

loop_macrosteps 5            # number of macro loops
loop_microsteps 40           # number of micro loops

tion1  Na
nion1  140                   # number of sodium atoms
tion2  Cl
nion2  100                   # number of chloride atoms
tion3  SO4
nion3  20

test_stable        yes
test_file          seawater.test
" > seawater.input

exe=./seawater
if [ -x $exe ]; then
  rm -fR state
  $exe
  rc=$?
  exit $rc
fi
exit 1

