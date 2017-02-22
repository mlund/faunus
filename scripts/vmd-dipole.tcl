# Script to draw dipole moment vectors for resid's
# and update for every frame in the tracjectory.
# Dipole moment based on the geometric center.
# Mikael Lund, September 2009.
#
# Usage:
# Specify resid range in for-loop below.
# $ source dipol.tcl
# $ enabletrace

proc vmd_draw_arrow {mol start end} {
  set end [vecadd $start [vecscale 0.4 $end]]
  set middle [vecadd $start [vecscale 0.8 [vecsub $end $start]]]
  set start [vecadd $start [vecscale -0.9 [vecsub $end $start]]]
  graphics $mol cylinder $start $middle radius 1.4
  #graphics $mol cylinder $start $end radius 1.4
  graphics $mol cone $middle $end radius 2.2
}

proc dipmom {args} {
  graphics 0 delete all
  for { set i 1 } { $i <= 1 } { incr i } {
    set sel1 [atomselect top "resid $i"]
    set muvec [measure dipole $sel1]
    set cm [measure center $sel1]
    #set vec [veclength $muvec]
    graphics 0 color red
    vmd_draw_arrow top $cm $muvec
  }
}

proc enabletrace {} {
  global vmd_frame
  #trace add variable ::vmd_frame([molinfo top]) write dipmom
  trace variable vmd_frame([molinfo top]) w dipmom
}

proc disabletrace {} {
  global vmd_frame
  #trace remove variable ::vmd_frame([molinfo top]) write dipmom
  trace vdelete vmd_frame([molinfo top]) w dipmom
}


