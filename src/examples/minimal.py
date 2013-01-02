from pyfaunus import *                   # import everything - prefer "import pyfaunus", though
cvar.atom.includefile("minimal.json")    # load atom properties
mcp = InputMap("minimal.input")          # open parameter file for user input
pot = Nonbonded_CoulombLJ_Cuboid(mcp)    # Define nonbonded pair-potential
spc = Space( pot.getGeometry() )         # create simulation space, particles etc.
salt = GroupAtomic(spc,mcp)              # group for salt particles
mv = AtomicTranslation(mcp,pot,spc)      # particle move class
mv.setGroup(salt)                        # tells move class to act on salt group
mv.move(100000)                          # move salt randomly 100000 times
print spc.info(), pot.info(), mv.info(), # final information
