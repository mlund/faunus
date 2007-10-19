#!/bin/bash
#
# This script will extract g(r) from the first two columns
# of a file, convert it into w(r) and adjust the offset so as
# to match that of a screened potential.
#
# M. Lund, October 2005

data=$1     # Data file
debye=30.44 # Debye length
beg=70      # Fitting range
end=140.     # -//-
QxQ=-5.87    # Charge product
lB=7.1      # Bjerrum length

function simple() {
echo "
set term table
plot '$data' u 1:($1*$1*($2-1))
" > tmp.gpi
gnuplot tmp.gpi | grep -v u
rm tmp.gpi
}

function virial() {
echo "
a=10.
b=4.
set term table
QxQ=39.6625
f(x)=7.13*QxQ/x*exp(-x/abs(a))+b
fit [$beg:$end] f(x) '$data' u 1:(-log(\$2)) via a,b
plot [0:$end] '$data' u 1:(-log(\$2)-b)
plot [$beg:600] f(x)-b
" > tmp.gpi
gnuplot tmp.gpi | grep -v u
#rm tmp.gpi                                                                               
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

virial

