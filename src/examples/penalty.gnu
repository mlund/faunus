set terminal svg enhanced fname "helvetica" fsize 20 size 1150,1050 dynamic
set output "penalty.svg"
set encoding utf8
set border lw 2
set multiplot
set size square 0.5,1
set lmargin 0
set rmargin 0
set view map
set xtics 0.5 border in scale 0,0 mirror norotate  offset character 0, 0.5, 0 autojustify
set ytics 0.5 border in scale 0,0 mirror norotate  offset character 0.4, 0, 0 autojustify
set ztics border in scale 0,0 nomirror norotate  offset character 0, 0, 0 autojustify
set rtics axis in scale 0,0 nomirror norotate  offset character 0, 0, 0 autojustify
set xrange [ -2.1 : 2.1 ] noreverse nowriteback
set yrange [ -2.1 : 2.1 ] noreverse nowriteback
set ylabel "y" offset 1,0
set xlabel "x" offset 0,1.04

set origin 0.02,0.26
unset colorbox
set cbrange [0:12]
set title "Initial histogram" font "helvetica, 25" offset 0, -0.7, 0
splot "histo1380" u 1:2:(log($3)) w p ps 1 pt 5 palette noti

set origin 0.47,0.26
set colorbox
set cbrange [0:12]
set cblabel "log(counts)"
set title "Final histogram" font "helvetica, 25" offset 0, -0.7, 0
splot "histo1681" u 1:2:(log($3)) w p ps 1 pt 5 palette noti

set origin 0.47,-0.24
set colorbox
set cbrange [0:18]
set cblabel "free energy"
set title "Penalty Function" font "helvetica, 25" offset 0, -0.7, 0
splot "pf_penalty" u 1:2:3 w p ps 1 pt 5 palette noti

set origin -0.0045,0.0015
set size 0.512,0.512
set yrange [0:18]
set ytics 4 nomirror 
set xtics offset 0,0.15
set xlabel "x" offset 0,0.7
set key spacing 2 at 2.1,6
set ylabel "free energy" offset 2,0 
set title "Final 1D profiles" font "helvetica, 25" offset 0, -0.7, 0
plot "histo1681" u 1:($2==0.5 ? log($3) : 1/0):(0.1) w circles lc rgb "#4992BA" fs transparent solid 0.60 border ti "log(histogram) at y=0.5",\
"pf_penalty" u 1:($2==0.5 ? $3 : 1/0):(0.1) w circles lc rgb "#00C795" fs transparent solid 0.60 border ti "penalty at y=0.5"

unset multiplot
