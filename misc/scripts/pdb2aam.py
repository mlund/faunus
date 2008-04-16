#!/sw/bin/python2.3

#
# Load atoms from PDB file and convert all
# residues to spheres, located at their center
# of mass (Amino Acid Model).
#
# NB: Uses the ScientificPython package
#
# M. Lund, 2008
#
# ...Husk: object.__dict__ !

from Scientific.IO.PDB import *
import sys

def atomWeight(atom):
  s=atom.__getitem__("name")
  #if atom.type=="ATOM":
  if s.find("O")==0: return 16.
  if s.find("C")==0: return 12.
  if s.find("N")==0: return 14.
  if s.find("S")==0: return 32.
  return 0.01;

# Calculate center-of-mass
def cm(res):
  p=[0,0,0,0]; wtot=0
  for atom in res:
    w=atomWeight(atom)
    wtot+=w
    for i in range(3):
      p[i]+=atom.position[i]*w
  for i in range(3): p[i]/=wtot
  p[3]=wtot
  return p

def findTitAtom(res):
  s="UNK"
  if res.name=="NTR":
    for atom in res:
      if atom.name=="N":
        return atom
  if res.name=="ASP": s="OD"
  if res.name=="GLU": s="OE"
  if res.name=="HIS": s="ND"
  if res.name=="TYR": s="OH"
  if res.name=="LYS": s="NZ"
  if res.name=="ARG": s="NH"
  if res.name=="CTR": s="OXT"
  for atom in res:
    if atom.name.find(s)==0:
      return atom
  return False

def aamString(res,atom):
  z=0.0
  mw=0
  r=2.0
  return res.name+" " + str(res.number) +" " \
      + str(atom.position[0])+" "+str(atom.position[1])+" "+str(atom.position[2])+" "\
      + str(z) + " " + str(mw) + " " + str(r)

# Load PDB file
model="MESO"
infile="test.pdb"
conf = Structure(infile)

# Loop over chains in PDB file
nchain=0
for chain in conf.peptide_chains:
  nchain+=1
  line=[]
  mw=0

  # Get amino acid spheres (if Mesoscopic model)
  if model=="MESO":
    for res in chain:
      oldname=res.name
      aa = cm(res)
      res.changeName("c"+oldname)
      line.append(
          res.name+" " + str(res.number) +" " \
        + str(aa[0])+" "+str(aa[1])+" "+str(aa[2])+" "\
        + "0.0" + " " + str(aa[3]) + " " + "3.5"  )
      res.changeName(oldname)

  # Search for titratable sites
  for res in chain:
    atom=findTitAtom(res)
    if atom!=False:
      line.append( aamString(res,atom) )

  # Take care of NTR and CTR
  first=0
  last=chain.__len__()-1
  chain[first].changeName("NTR")
  chain[last].changeName("CTR")
  for n in [first, last]:
    atom=findTitAtom(chain[n])
    if atom!=False:
      line.append( aamString(chain[n],atom) )

  # Print AAM file
  f=open(infile+"-"+str(nchain)+".aam","w")
  f.write( "# "+infile+" - Chain "+str(nchain) \
      + " ("+str(chain.__len__())+" residues)\n" \
      + str(line.__len__())+"\n" )
  for l in line:
    f.write(l+"\n")
  f.close()

sys.exit(0)

