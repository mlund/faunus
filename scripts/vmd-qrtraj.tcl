# Used for visualizing charge and radius changes in VMD trajectories.
# Loads a charge-radius text files where each line corresponds
# to one frame with alternating charge - white space - radius values for each
# atom.
#
# Inspired by http://klein-group.icms.temple.edu/cpmd-vmd/files/ar3plus-charge.vmd 
#
# Usage:
#
# (1) Load trajectory into VMD - number of frames must match 'qrtraj.dat' file
# (2) in the VMD console, `source vmd-qrtraj.tcl`

set molid 0
set n [molinfo $molid get numframes]
set fp [open "qrtraj.dat" r]
for {set i 0} {$i < $n} {incr i} {
    set qrdata($i) [gets $fp]
}
close $fp

proc do_update {args} {
    global qrdata molid
    set a [molinfo $molid get numatoms]
    set f [molinfo $molid get frame]
    for {set i 0} {$i < $a} {incr i 1} {
        set s [atomselect $molid "index $i"]
        $s set charge [lindex $qrdata($f) [expr 2*$i] ]
        $s set radius [lindex $qrdata($f) [expr 2*$i+1 ] ]
  }
}

trace variable vmd_frame($molid) w do_update
mol colupdate 0 $molid on
animate goto start
do_update
