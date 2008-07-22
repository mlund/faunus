#!/bin/bash
#
# This script will extract g(r) from the first two columns
# of a file, convert it into w(r) and adjust the offset so as
# to match that of a screened potential.
#
# M. Lund, October 2005

data=$1     # Data file
debye=13.4  # Debye length
beg=40      # Fitting range
end=60.     # -//-
QxQ=-10.    # Charge product
lB=7.1      # Bjerrum length

function simple() {
echo "
set term table
plot '$data' u 1:(-log(\$2))
" > tmp.gpi
gnuplot tmp.gpi | grep -v u
rm tmp.gpi
}

function debyehuckel() {
echo "
set term table
f(x)=$lB*$QxQ/x*exp(-x/$debye)+a
fit [$beg:$end] f(x) '$data' u 1:(-log(\$2)) via a
plot '$data' u 1:(-log(\$2)-a), f(x)-a
" > tmp.gpi
gnuplot tmp.gpi | grep -v u
rm tmp.gpi
}

debyehuckel

