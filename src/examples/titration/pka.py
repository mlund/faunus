#!/usr/bin/env python
import faunus

print faunus.faunus_splash(),            # Faunus spam
inp  = faunus.inputfile("pka.conf")      # read input file
loop = faunus.mcloop(inp)                # set Markov chain loop lengths
nvt  = faunus.canonical()                # use the canonical ensemble
con  = faunus.cell(inp)                  # use a spherical simulation container...
pot  = faunus.interaction_hscoulomb(inp) # ... and a Coulomb + hard sphere potential
sm   = faunus.saltmove(nvt,con,pot,inp)  # create object for salt movements
prot = faunus.macromolecule()
aam  = faunus.ioaam(con.atom)            # load/save configurations from/to disk
pqr  = faunus.iopqr(con.atom)            # load/save configurations from/to disk
prot.add(con, inp)                                # laod protein
prot.move(con, -prot.cm);                # ..translate it to origo (0,0,0)
prot.accept(con)                         # ..accept translation
salt = faunus.salt()                     # Group for salt and counter ions
salt.add(con, inp)                       #   Insert sodium ions
aam.load(con, "confout.aam")             # Load old config (if present)

tit = faunus.chargereg(nvt,con,pot,salt,inp.getflt("pH",7.))

sys = faunus.systemenergy(pot.energy(con.p)) # System energy analysis
print con.info(), tit.info(), pot.info(), con.atom.info(),

while(loop.macroCnt()):                  # Markov chain 
    while (loop.microCnt() ):
        sys+=sm.move(salt)               # Displace salt particles

        sys+=tit.titrateall()            # Titrate protein sites
        prot.charge(con.p)               # Re-calc. protein charge
        prot.dipole(con.p)               # Re-calc. dipole moment
    sys.update(pot.energy(con.p))        # Update system energy
    aam.save("confout.aam", con.p)       # Save config. to disk
    pqr.save("confout.pqr", con.p, tit)  # Save PQR file to disk - cool in VMD!
    print loop.timing(),                 # Show progress

print sys.info(), sm.info(), tit.info(), salt.info(con),
print prot.info(), loop.info(),          # Print final results

