#!/usr/bin/env python

import faunus

# show Faunus information
print faunus.faunus_splash()

# read input file
inp = faunus.inputfile('widom.conf')

# set Markov chain loop lengths
loop = faunus.mcloop(inp)

# use the canonical ensemble
nvt = faunus.canonical()

# use a spherical simulation container...
cell = faunus.cell(inp)

# ... and a Coulomb + hard sphere potential
pot = faunus.interaction_hscoulomb(inp)

# create object for salt movements
sm = faunus.saltmove(nvt, cell, pot, inp)

# define some groups for mobile ions
salt = faunus.salt()

# insert some ions
salt.add(cell,inp)

# file I/O object, read initial config from disk (if present)
aam = faunus.ioaam(cell.atom)
aam.load(cell, 'widom.aam')

# virial analysis
# currently missing

# object for multiple particle insertion
wid1 = faunus.widom(10)

# object for single particle insertion with charge scaling
#wid2 = faunus.widomSW(10)

# detect all species in the cell for both
wid1.add(cell)
#wid2.add(cell)

# track system energy
sys = faunus.systemenergy(pot.energy(cell.p))

# print initial information
print cell.info()
print cell.atom.info()
print pot.info()
print salt.info()
print inp.info()

# Markov chain
while(loop.macroCnt()):
    while(loop.microCnt()):
        sys += sm.move(salt)        # displace salt particles
        wid1.insert(cell, pot)      # Widom particle insertion analysis
        #wid2.insert(cell, pot)      #  - // -
        #virial.sample(cell, pot)    # virial sampling
    sys.update(pot.energy(cell.p))  # update system energy
    aam.save('widom.aam', cell.p)   # save particle configuration to disk
    print loop.timing()             # show progress
    
# print final information and results
print sys.info()
print sm.info()
print wid1.info()
print loop.info()
