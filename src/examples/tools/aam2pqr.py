#!/usr/bin/env python

import sys
import os.path
import faunus as fau

if (len(sys.argv)==3):

  infile  = sys.argv[1]
  outfile = sys.argv[2]
  if (os.path.isfile(infile)==False):
    print "Input structure does not exist."
    sys.exit(0)

  print fau.faunus_splash()
  con  = fau.cell(1000)
  prot = fau.group()
  aam  = fau.ioaam()
  pqr  = fau.iopqr()

  prot.add(con, aam.load(infile))
  pqr.save(outfile, con.p)

else:

  print "Usage: aam2pqr infile.aam outfile.pqr"
