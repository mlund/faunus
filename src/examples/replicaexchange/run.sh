#!/bin/bash

function mkinput() {
echo "boxlen        $box
temperature   298.15
bjerrum       $bjerrum
LJeps         0.30
nprot1        5
protein1      test.mol1.aam
nprot2        5
protein2      test.mol2.aam
mrdp          2
mtdp          100
ctdp          10
dppercent     0.0008 
e_r           80
debyelen      100
moltrans_dp   $tdp
isobar_pressure      $pressure
isobar_dp            $voldp
isobar_maxlen        500
isobar_minlen        50
isobar_binlen        1
isobar_penalize      no
isobar_penalty       0.5
isobar_scalepen      1 
penaltyupdate        no
atomfile      faunatoms.dat
volr          1
tr            1
rr            1
clt           1 
tr_rot        yes
movie         no
" > $id.conf
}

function mkrepinput {
echo "
macrosteps $macro
microsteps $micro
temper     $temper
" > replica.conf
}

macro=10
pressure=0.000001

# 1. Generate replica input files
for id in a b c
do
  if [ "$id" == "a" ]; then box=150; bjerrum=12;  tdp=10; voldp=0.7; fi
  if [ "$id" == "b" ]; then box=200; bjerrum=10;  tdp=100; voldp=1; fi
  if [ "$id" == "c" ]; then box=250; bjerrum=4;   tdp=200; voldp=2; fi
  mkinput
done

export OMP_NUM_THREADS=1

# 2. Eq run
rm -f *.dump
macro=10
micro=10000
temper="yes"
mkrepinput
./replica
exit

# 3. Production
micro=10000
temper="yes"
mkrepinput
./replica

