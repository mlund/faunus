# Script to read residue charges per frame from file
# Chris Evers, March 2011
#
# Inspired by http://klein-group.icms.temple.edu/cpmd-vmd/files/ar3plus-charge.vmd 
#
# Usage:
# $ source vmd-qtraj.tcl

set molid 0
# read in charge data
set n [molinfo $molid get numframes]
puts "reading charges"
set fp [open "q.traj" r]
for {set i 0} {$i < $n} {incr i} {
  set chrg($i) [gets $fp]
}
close $fp

# procedure to change the charge field from the data in $chrg
proc do_charge {args} {
  global chrg molid
  set a [molinfo $molid get numatoms]
  set f [molinfo $molid get frame]
  for {set i 0} {$i < $a} {incr i} {
    set s [atomselect $molid "index $i"]
    $s set charge [lindex $chrg($f) $i]
  }
}

# turn on update of the charge info and coloring in each frame.
trace variable vmd_frame($molid) w do_charge
mol colupdate 0 $molid on

# set color mapping parameters
mol scaleminmax 0 $molid -1.0 1.0
color scale method RGB
color scale midpoint 0.25
color scale min 0.0
color scale max 1.0
# go back to the beginning and activate the additional features.
animate goto start
do_charge
