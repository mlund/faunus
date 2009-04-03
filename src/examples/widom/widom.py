#!/usr/bin/env python
import faunus

print faunus.faunus_splash(),            # show Faunus information
inp  = faunus.inputfile('widom.conf')    # read input file
loop = faunus.mcloop(inp)                # set Markov chain loop lengths
nvt  = faunus.canonical()                # use the canonical ensemble
cell = faunus.cell(inp)                  # use a spherical simulation container...
pot  = faunus.interaction_hscoulomb(inp) # ... and a Coulomb + hard sphere potential
sm   = faunus.saltmove(nvt,cell,pot,inp) # create object for salt movements
wid  = faunus.widom(10)                  # object for multiple particle insertion
salt = faunus.salt()                     # define some groups for mobile ions
salt.add(cell,inp)                       #    insert some ions
aam  = faunus.ioaam(cell.atom)           # load/save configurations from/to disk
aam.load(cell, 'widom.aam')              #    get initial configuration (if present)
wid.add(cell)                            # detect all species in the cell
sys  = faunus.systemenergy(pot.energy(cell.p)) # track system energy

# print initial information
print cell.info(), pot.info(), salt.info(), inp.info()

# begin simulation
while(loop.macroCnt()):
    while(loop.microCnt()):
        sys += sm.move(salt)        # displace salt particles
        wid.insert(cell, pot)       # Widom particle insertion analysis
    sys.update(pot.energy(cell.p))  # update system energy
    aam.save('widom.aam', cell.p)   # save particle configuration to disk
    print loop.timing(),            # show progress
    
# print final information and results
print sys.info(), sm.info(), wid.info(), loop.info()
