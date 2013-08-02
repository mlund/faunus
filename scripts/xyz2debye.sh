#!/bin/bash

# This will take an XYZ trajectory file, split it and pass each
# frame on to "Debyer", https://code.google.com/p/debyer to
# calculate the isotropic structure factor using the Debye formula.
#
# The expected xyz format is:
#
# line 1:   number of particles, N
# line 2:   molname
# line 3:   atom x y z (coordinates in nanometers!)
# line N+2: next xyz frame
#
# The "atom" type is ignored in this script as we calculate
# the structure factor and particle form factors are set to
# unity.

# User input
in=$1      # .xyx trajectory file as 1st argument
box=0      # cube side length (nm)
qmax=5.0   # maximum q value (nm)
dq=0.01    # q spacing (nm)

read -p "Box length (nm): " box                                       

debyer=$HOME/Downloads/debyer-r77/debyer/debyer

# Calculated variables
cutoff=`python -c "print $box/2.01"`
N=`head -n 1 $in`
qmin=`python -c "print 2*3.1416/($box/2.)"`

echo "XYZ trajectory file     "$in
echo "Number of particles     "$N
echo "Cube side length (nm)   "$box
echo "Cutoff distance (nm)    "$cutoff
echo "q range (1/nm)          ["$qmin":"$qmax"] dq="$dq
echo "Debyer executable       "$debyer

echo 
read -p "Continue? "

prefix="_cm."
rm -fR _cm.*

# First split xyz trajectory to individual files
split -a 4 -l `python -c "print $N+2"` $in $prefix

# loop over each file and calc. S(q)
files=`ls $prefix*`
for f in $files
do
  mv $f tmp.xyz
  $debyer --sf -a $box -b $box -c $box -r $cutoff -f $qmin -t $qmax -s $dq -o $f.dat tmp.xyz 
  rm tmp.xyz
done

# average all frames to final S(q)
LANG=en
awk '{ sum[$1]+=$2; cnt[$1]++ } END { for (i in sum) print i, sum[i]/cnt[i] | "sort -n" }' $prefix*.dat > sofq.dat

