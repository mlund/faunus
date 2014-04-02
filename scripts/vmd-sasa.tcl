#
# VMD script to measure SASA for residues in proteins
#
# This is a VMD script to measure solvent accessible
# surface areas (SASA) of residues in proteins. The
# individual SASAs (in angstrom squared) is saved
# to a one column file, `out.sasa`, where each row
# corresponds to a residue. The probe radius is 1.5
# angstrom.
#
# Usage:
#
#    $ vmd -dispdev text infile.pdb -e vmd-sasa.tcl
#
set sum 0
set outfile [open "out.sasa" w]
set sel [atomselect top "protein"]
set residues [lsort -integer -unique [$sel get resid]]
foreach r $residues {
  set area [measure sasa 1.5 $sel -restrict [atomselect top "resid $r"]]
  set sum [expr "$sum + $area"]
  puts $outfile $area
}
puts "Total SASA = $sum"
close $outfile
quit
