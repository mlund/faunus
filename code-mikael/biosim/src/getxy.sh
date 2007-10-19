#!/bin/bash
#
# M. Lund, October 2005

data=$2     # Data file
xy=$1       # Column to extract

function plot() {
echo "
set term table
plot '$data' u ${xy}
" > tmp.gpi
gnuplot tmp.gpi | grep -v u | grep -v nan
rm tmp.gpi
}

plot
