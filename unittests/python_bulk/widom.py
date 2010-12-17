#!/usr/bin/env python
import pyfaunus as fau

print fau.faunus.splash(),                # faunus spam
inp  = fau.inputfile('widom.conf')        # read input file
loop = fau.mcloop(inp)                    # set Markov chain loop lengths
nvt  = fau.canonical()                    # use the canonical ensemble
con  = fau.box(inp)                       # use a spherical simulation container...
pot  = fau.interaction_hsminimage(inp)    # ... and a Coulomb + hard sphere potential
sm   = fau.saltmove(nvt,con,pot,inp)      # create object for salt movements
wid  = fau.widom(10)                      # object for multiple particle insertion
widSW= fau.widomSW(10)                    # object for single particle insertion w. charge scaling (HS-Coulomb only!)
salt = fau.salt()                         # define some groups for mobile ions
salt.add(con,inp)                         #    insert some ions
aam  = fau.ioaam()                        # load/save configurations from/to disk
aam.load(con, 'widom.aam')                #    get initial configuration (if present)
wid.add(con)                              # detect all species in the container
widSW.add(con);                           # - // -
sys = fau.systemenergy(pot.energy(con.p)) # track system energy

print inp.info() + con.info() + pot.info() + salt.info()

while(loop.macroCnt()):                   # markov chain
    while(loop.microCnt()):
        sys += sm.move(salt)              # displace salt particles
        wid.insert(con, pot)              # Widom particle insertion analysis
        widSW.insert(con, pot)            # - // -
    sys.update(pot.energy(con.p))         # update system energy
    aam.save('widom.aam', con.p)          # save particle configuration to disk
    print loop.timing(),                  # show progress
    
print sys.info() + sm.info() + wid.info() + widSW.info() + loop.info()
