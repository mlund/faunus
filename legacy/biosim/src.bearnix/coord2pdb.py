#!/sw/bin/python2.3

from Scientific.IO.ArrayIO import *
from Scientific.IO.PDB import *

d = readFloatArray(".coord.a")
f = PDBFile("test.pdb", "w")

for i in range(len(d)):
  f.nextResidue("UNK")
  f.writeAtom("O", d[i])

f.close()
