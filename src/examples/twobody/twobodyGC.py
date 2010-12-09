#!/usr/bin/env python

import sys, random
import pyfaunus as fau
from math import *

print fau.faunus_splash()                # Faunus info
shout = fau.notifyUser()                 # A class for notifying user about progress
inp   = fau.inputfile("twobody.conf")    # Read input file
test  = fau.checkValue(inp)              # Test output
loop  = fau.mcloop(inp)                  # Set Markov chain loop lengths
cell  = fau.cell(inp)                    # We want a spherical, hard cell
nmt   = fau.grandcanonical()
pot   = fau.interaction_coulomb(inp)     # Functions for interactions
dst   = fau.distributions()              # Distance dep. averages
gofr_pp = fau.atomicRdf(0.5,200)         # Amino acid rdf between the two proteins

# IO
io  = fau.io()
pqr = fau.iopqr()                        # PQR output (pos, charge, radius)
xtc = fau.ioxtc(1000.)                   # Gromacs xtc output
aam = fau.ioaam()                        # Protein input file format is AAM

gofr_pp_rf = inp.getflt("gofr_pp_rf",0)  # How often should gofr_pp be sampled?
g = fau.vector_macromolecule()           # Group for proteins

# Protein setup
mr = fau.macrorot(nmt, cell, pot)        # Class for macromolecule rotation
dm = fau.dualmove(nmt, cell, pot)        # Class for 1D macromolecular translation
dm.setup(inp)
dm.load( inp, g, 80.)                    # Load proteins and separate them 
vg = fau.vector_group()                  # Make a vector of just the proteins
vg.push_back(g[0])                       # (used for trajectory output)
vg.push_back(g[1])

# Salt setup
salt = fau.salt()                        # Group for mobile ions
salt.add(cell, inp)                      #   Add salt particles
sm = fau.saltmove(nmt, cell, pot, inp)   # Class for salt movements

# Titration and GC moves
if (inp.getflt("tit_runfrac", 0.5)<1e-5 ):   # Should we disable titration?
  for i in range( len(fau.cvar.atom.list) ):
    atom.list[i].pka=0

sb  = fau.saltbath(nmt,cell,pot,inp,salt)    # Class for GC salt insertion
tit = fau.GCchargereg(nmt,cell,pot,inp)      # Class for proton titration

if (nmt.load(cell, "gcgroup.conf")==True):                
  aam.load(cell,"confout.aam")    # Read initial config. from disk (if present)
  g[0].masscenter(cell)           # Load old config (if present)
  g[1].masscenter(cell)           # ...and recalc mass centers

sysu = fau.systemenergy(          # Track total system energy
    pot.energy(cell.p, g[0], g[1]) +
    pot.energy(cell.p, salt) +
    pot.internalElectrostatic(cell.p, g[0]) +
    pot.internalElectrostatic(cell.p, g[1]) +
    pot.internal(cell.p, salt) )

print( "# ---- Initial information ----" +
    inp.info() + cell.info() + tit.info() + pot.info() )

print "# ---- Runtime output ----"
while loop.macroCnt():                    # Markov chain 
  while loop.microCnt():
    rnd = random.randint(1,4)
    if (rnd==1):
      if (random.random()>0.5):
        sysu += sm.move(salt)                # Translate salt
      else:
        for i in range ( salt.size() ):
          sysu+=sb.move()                    # Grand canionical salt
    if (rnd==2):
      for i in range(2):
        j = random.randint(0,1)
        sysu += mr.move(g[j])                # Rotate proteins
    if (rnd==3):
      sysu += dm.move(g[0], g[1])            # Translate both proteins
    if (rnd==4):
      sysu += tit.titrateall()               # Titrate sites on the protein
      dst.add( "Q1",dm.r, g[0].charge(cell.p) ) 
      dst.add( "Q2",dm.r, g[1].charge(cell.p) ) 
      dst.add( "MU1",dm.r,g[0].dipole(cell.p) ) 
      dst.add( "MU2",dm.r,g[1].dipole(cell.p) ) 

    if ( random.random()>0.9 ):
      g[0].dipole(cell.p)
      g[1].dipole(cell.p)
      dst.add( "protein1 z-dipcomp", dm.r, g[0].mu.z/g[0].mu.len() )
      dst.add( "protein2 z-dipcomp", dm.r, g[1].mu.z/g[1].mu.len() )

    if ( random.random()>0.95 ):                    # Track system energy
      dst.add("Utot", dm.r, pot.energy(cell.p)) 

    if ( random.random()<gofr_pp_rf ):              # Track protein-protein g(r)
      gofr_pp.update( cell.p, g[0], g[1] ) 

    if ( random.random()>.99 ):                     # Save trajectory of proteins
      xtc.save("coord.xtc", cell.p, vg) 

    # End of inner loop
    
  sysu.update(                                  # Update system energy
    pot.energy( cell.p, g[0], g[1] ) +
    pot.energy( cell.p, salt ) +
    pot.internalElectrostatic( cell.p, g[0] ) +
    pot.internalElectrostatic( cell.p, g[1] ) +
    pot.internal( cell.p, salt) )               # System energy analysis

  gofr_pp.write("rdfatomic.dat")              # Write interprotein g(r) - Atomic
  dm.gofr.write("rdfprot.dat")                # Write interprotein g(r) - CM
  dst.write("distributions.dat")              # Write other distributions
  aam.save("confout.aam", cell.p)             # Save config. for next run

  tit.applycharges(cell.trial)                # Set average charges on all titratable sites
  pqr.save("confout.pqr", cell.trial)         # ... save PQR file
  pqr.save("confout_ns.pqr", cell.trial,vg)   # ... save PQR file (proteis only)
  aam.save("confout_smeared.aam", cell.trial) 
  cell.trial=cell.p                           # Restore original charges!

  io.writefile("gcgroup.conf", nmt._print()) 

  cell.check_vector()                         # Check sanity of particle vector
  print loop.timing(),                        # Show progress

  # End of outer loop

xtc.close()                                   # Close xtc file for writing

print( "# ---- Final information ----"
     + loop.info() + salt.info(cell)
     + sm.info() + mr.info() + dm.info()
     + sysu.info() + g[0].info() + g[1].info() + cell.info()
     + tit.info() + sb.info() )

# Unit tests
sysu.check(test) 
sb.check(test) 
dm.check(test) 
sm.check(test) 
mr.check(test) 
test.check("Protein1_charge", g[0].Q.avg()) 
test.check("Protein2_charge", g[1].Q.avg())
print test.report() 

shout.message("TwobodyGC finished!",False) 
sys.exit( test.returnCode() )

