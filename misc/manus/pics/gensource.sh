#!/bin/bash
enscript --color -C -o s.ps -E -B -2r --font=CourierBold6 -M Env10 -L 26 $1
ps2pdf14 s.ps
pdfcrop s.pdf


