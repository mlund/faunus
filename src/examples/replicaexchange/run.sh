#!/bin/bash

function mkinput() {
echo "boxlen        100
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
debyelen      40

moltrans_dp   100

isobar_pressure      0.000001
isobar_dp            0.40
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

for id in a b c
do
  rm -f *.dump
  if [ "$id" == "a" ]; then bjerrum=7; fi
  if [ "$id" == "b" ]; then bjerrum=4; fi
  if [ "$id" == "c" ]; then bjerrum=1; fi
  mkinput
done

./replicaex

