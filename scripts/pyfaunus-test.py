# DYLD_LIBRARY_PATH="/Users/mikael/miniconda/lib/" python test.py

import pyfaunus as fau
import numpy as np

mcp = fau.InputMap('minimal.json')
spc = fau.Space(mcp)

g = spc.groupList()[0]
cm = fau.massCenter(spc.geo(), spc.p(), g)
print cm.x, cm.y, cm.z
a = np.array( cm )
print a