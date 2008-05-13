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

# ---------- FUNCTIONS: -----------

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
    if res[0].name=="N":
      return res[0]
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

def aamStr(resname,number,pos,z,mw,r):
  return resname+" " + str(number) +" " \
      + str(pos[0])+" "+str(pos[1])+" "+str(pos[2])+" "\
      + str(z) + " " + str(mw) + " " + str(r)

def aamString(res,atom):
  z=0.0
  mw=0
  r=2.0
  return res.name+" " + str(res.number) +" " \
      + str(atom.position[0])+" "+str(atom.position[1])+" "+str(atom.position[2])+" "\
      + str(z) + " " + str(mw) + " " + str(r)

# get Name, charge and weight
def findHetatom(name):
  r=False
  if name=="CA": r=[2,40.078]
  if name.find("V")==0: r=[4,50.94]   # III,IV,V
  #if name.find("MG")==0: r=[2,24.3]
  if name.find("FE")==0: r=[2,55.85]  # II,III
  if name.find("CO")==0: r=[3,58.93]
  if name.find("CU")==0: r=[2,63.55]  # I,II
  if name.find("CR")==0: r=[3,51.996]
  if name.find("ZN")==0: r=[2,65.409]
  if name.find("MN")==0: r=[4,54.938] # II,IV
  return r

# ---------- MAIN PROGRAM -----------

# Load PDB file
model="DUMMY"
infile="1jb0.pdb"
conf = Structure(infile)

# Loop over chains in PDB file
nchain=0  # Chain counter
aam=[]    # Array of aam lines
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
  last=len(chain)-1
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

  # Add AAM info to array
  aam.append(line[:])

# Write an AAM file with ALL chains and relevant HETATOM's
# Search for HETATOM's
line=[]
for res in conf:
  for atom in res:
    if atom.type()=="HETATM":
      name=atom.__getitem__("name")
      het=findHetatom(name)
      if het!=False:
        line.append( aamStr(name,res.number,atom.position,het[0],het[1],3.0) )
aam.append(line[:])
n=0
for chain in aam:
  n+=len(chain)
print n
for chain in aam:
  for line in chain:
    print line
