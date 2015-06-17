import numpy as np
import json, sys, os
from subprocess import call, check_output
from shutil import copyfile

pfx = os.path.join( os.path.dirname(__file__), "polymers")
try:
  copyfile( pfx+'.test', 'polymers.test' )
  copyfile( pfx+'.state', 'state' )
except: pass

d = {
    "atomlist" : {
      "Na" : { "q": 1.0, "r":2.0, "dp":110 },
      "Cl" : { "q":-1.0, "r":2.0, "dp":110 },
      "MM" : { "q": 1.0, "r":3.0, "dp":10 }
      },

    "moleculelist" : {
      "counterions" : { "atoms":"Cl", "atomic":True, "Ninit":32 },
      "polymer" : {
        "structure":"polymers.aam", "Ninit":8,
        "bonds" : {
          "0 1" : { "type":"harmonic", "k":0.0557, "req":0.0 },
          "1 2" : { "type":"harmonic", "k":0.0557, "req":0.0 },
          "2 3" : { "type":"harmonic", "k":0.0557, "req":0.0 }
          }
        }
      },

    "energy" : {
      "nonbonded" : { "coulomb" : { "epsr" : 78.7 } }
      },

    "moves" : {
      "atomtranslate" : {
        "counterions" : { "prob":1.0, "peratom":True },
        "polymer"     : { "prob":1.0, "permol":True, "peratom":True }
        },
      "moltransrot" : {
        "polymer" : { "dp":50, "dprot":6, "permol":True, "dir":"1 1 1" } 
        },
      "crankshaft" : {
        "polymer" : { "dp":6.3, "minlen":2, "maxlen":10, "permol":True }
        },
      "pivot" : {
        "polymer" : { "dp":6.3, "minlen":1, "maxlen":1000, "permol":True }
        },
      "isobaric" : { "dp":2, "pressure":113.2, "prob":0.1 }
      },

    "system" : {
      "cuboid"       : { "len" : 200 },
      "sphere"       : { "radius" : 100 },
      "mcloop"       : { "macro":10, "micro":20000 },
      "unittest"     : { "testfile":"polymers.test", "stable":False }
      }
    }

print >> open('polymers.json', 'w+'), json.dumps(d, indent=4)

# molecule with four atoms
print """4
MM  0   0.00   0.00   0.00    1.0   1  3.0
MM  1   7.60   0.00   0.00    1.0   1  3.0
MM  2   0.00   7.60   0.00    1.0   1  3.0
MM  3   7.60   7.60   0.00    1.0   1  3.0"""

