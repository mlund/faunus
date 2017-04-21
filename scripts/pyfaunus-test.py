# DYLD_LIBRARY_PATH="/Users/mikael/miniconda/lib/" python test.py

import pyfaunus as mc
import numpy as np

mcp    = mc.InputMap('../src/examples/minimal.json')
spc    = mc.Space(mcp)
geo    = spc.geo
groups = spc.groupList()

print "number of particles            = ", len(spc.p)
print "number of groups               = ", len(groups)
print "system volume                  = ", geo.getVolume()

cm = mc.massCenter(spc.geo, spc.p, groups[0])
print "mass center                    = ", np.array(cm)
print "distance between two particles = ", geo.dist(spc.p[0], spc.p[1])

print "index in 1st group:\n", groups[0].range()