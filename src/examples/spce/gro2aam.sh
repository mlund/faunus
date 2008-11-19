#!/bin/bash

#Gro:  res name atomnr x y z
#AAM:  name atomnr x y z q mw r

#genbox -cs spc216 -o water.gro -box $len $len $len

cat "$1" | sed s/HW1/"HW1 "/ | sed s/HW2/"HW2 "/ | sed s/OW/"OW "/ | gawk '{print $2" "$3" "$4*10-50" "$5*10-50" "$6*10-50" 0 1 1.0"}'

