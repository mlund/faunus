#!/usr/bin/env python

import sys
import pyfaunus as faunus

print faunus.faunus_splash(),            # Faunus spam
config = "pka.conf"                      # Default input (parameter) file
if len(sys.argv) == 2:
    config = sys.argv[1]                 # ..also try to get it from the command line
inp  = faunus.inputfile(config)          # read input file
loop = faunus.mcloop(inp)                # set Markov chain loop lengths
nvt  = faunus.canonical()                # use the canonical ensemble
con  = faunus.cell(inp)                  # use a spherical simulation container...
pot  = faunus.interaction_hscoulomb(inp) # ... and a Coulomb + hard sphere potential
sm   = faunus.saltmove(nvt,con,pot,inp)  # create object for salt movements
prot = faunus.macromolecule()
aam  = faunus.ioaam()                    # load/save configurations from/to disk
pqr  = faunus.iopqr()                    # load/save configurations from/to disk
prot.add( con,
        aam.load(inp.getstr("protein"))) # load protein structure
prot.masscenter(con)
prot.move(con, -prot.cm)                 # ..translate it to origo (0,0,0)
prot.accept(con)                         # ..accept translation
salt = faunus.salt()                     # Group for salt and counter ions
salt.add(con, inp)                       #   Insert sodium ions
aam.load(con, "confout.aam")             # Load old config (if present)

tit = faunus.chargereg(nvt, con, pot, salt, inp.getflt("pH",7.))

sys = faunus.systemenergy(pot.energy(con.p)) # System energy analysis
print inp.info(), con.info(), tit.info(), pot.info(),

while loop.macroCnt():                   # Markov chain 
    while loop.microCnt():
        sys += sm.move(salt)             # Displace salt particles
        sys += tit.titrateall()          # Titrate protein sites
        prot.charge(con.p)               # Re-calc. protein charge
        prot.dipole(con.p)               # Re-calc. dipole moment
    sys.update(pot.energy(con.p))        # Update system energy
    aam.save("confout.aam", con.p)       # Save config. to disk
    pqr.save("confout.pqr", con.p, tit)  # Save PQR file to disk - cool in VMD!
    print loop.timing(),                 # Show progress

print sys.info(), sm.info(), tit.info(), salt.info(con),
print prot.info(), loop.info(),          # Print final results

