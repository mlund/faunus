#!/sw/bin/python2.3

#
# Load atoms from PDB file and convert all
# residues to spheres, located at their center
# of mass (Amino Acid Model).
#
# NB: Uses the ScientificPython package
#
# M. Lund, 2005
#
# ...Husk: object.__dict__ !

from Scientific.IO.PDB import *
import sys

def charge(name):
  return 0.

# Approximate radius from weight
def radius(mw):
  rho=0.9
  return pow( mw/rho * 3./4./3.14*1.661, 0.3333 )

# Return weight of an atom
def weight(atomname):
  if atomname.find("O")==0: return 16.
  if atomname.find("C")==0: return 12.
  if atomname.find("N")==0: return 14.
  if atomname.find("S")==0: return 32.
  return 1.0

# Calculate center-of-mass
def cm(residue):
  p=[0,0,0,0]; wtot=0
  for atom in residue:
    w=weight(atom.name)
    wtot+=w
    for i in range(3): p[i]+=atom.position[i]*w
  for i in range(3): p[i]/=wtot
  p[3]=wtot
  return p

conf = Structure('test.pdb')
print conf
cnt=0
for res in conf.residues:
  cnt+=1;
  if (cnt-res.number)>5:
    print cnt, res.number;
    cnt=1

sys.exit(0)
for residue in conf.residues:
  if residue.name!="HOH":
    aa = cm(residue)
    #print residue.name, residue.number, aa[0], aa[1], aa[2], aa[3], radius(aa[3])
    print "%4s %5d %6.2f %6.2f %6.2f %4.1f %4.0f %4.1f" % \
        (residue.name, residue.number, aa[0], aa[1], aa[2], charge(residue.name), aa[3], radius(aa[3]))

