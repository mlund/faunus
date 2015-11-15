# DYLD_LIBRARY_PATH="/Users/mikael/miniconda/lib/" python test.py

import pyfaunus as fau
import numpy as np

mcp = fau.InputMap('minimal.json')
spc = fau.Space(mcp)
print spc.info()

for i in spc.p():
  print i.charge

