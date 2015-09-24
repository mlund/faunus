import numpy as np
import json, sys, os
from subprocess import call, check_output
from shutil import copyfile

pfx = os.path.join( os.path.dirname(__file__), "polyelectrolytes")
#try:
#  copyfile( pfx+'.state', 'state' )
#  copyfile( pfx+'.test', 'polymers.test' )
#except: pass

d = {
    "atomlist" : {
      "Na" : { "q": 1.0, "r":2.0, "dp":110, "activity":0.001 },
      "Cl" : { "q":-1.0, "r":2.0, "dp":110, "activity":0.001 },
      "HA" : { "q": 0.0, "r":2.0, "dp":4.0},
      "A-" : { "q":-1.0, "r":2.0, "dp":4.0},
      "ALA": { "q": 1.0, "r":2.0, "dp":4.0 }
      },
    "processes": {
        "HA-A": { "pKd": 4.8,   "pX": -7.8, "bound": "HA", "free": "A-" }
    },

    "moleculelist" : {
      "salt" : { "atoms":" Na Cl", "atomic":True, "Ninit":10 },
      "polymer" : {
        "Ninit"        :10, 
        "blockpolymer":"HA",
        "atoms"       :"HA",
        "blocksizes"  :"2",
        "req"         :"6.0",
        "k"           :"1.0"
#         "fasta" : "AAAA", "Ninit":8
#        "structure":"polymers.aam", "Ninit":8,
#        "bonds" : {
#          "0 1" : { "type":"harmonic", "k":0.0557, "req":4.0 },
#          "1 2" : { "type":"harmonic", "k":0.0557, "req":4.0 },
#          "2 3" : { "type":"harmonic", "k":0.0557, "req":4.0 }
#          }
        }
      },

    "energy" : {
      "eqstate"   : { "processfile": "polyelectrolyte.json" },
      "nonbonded" : { "coulomb" : { "epsr" : 78.7 } }
      },

    "moves" : {
      "atomtranslate" : {
        "salt" : { "prob":0.1, "peratom":True },
        "polymer"     : { "prob":1.0, "permol":True, "peratom":True }
        },
      "gctit"    : { "molecule": "salt", "prob": 1.00 
        },
      "moltransrot" : {
        "polymer" : { "prob":0.1, "dp":50, "dprot":6, "permol":True, "dir":"1 1 1" } 
        },
#      "crankshaft" : {
#        "polymer" : { "dp":6.3, "minlen":2, "maxlen":10, "permol":True, "prob":1.0 }
#        },
#      "pivot" : {
#        "polymer" : { "dp":6.3, "minlen":1, "maxlen":10, "permol":True, "prob":1.0 }
#        },
      "isobaric" : { "dp":0.5, "pressure":2.3, "prob":0.5 }
      },

    "system" : {
      "cuboid"       : { "len" : 2000 },
      "sphere"       : { "radius" : 100 },
      "mcloop"       : { "macro":10, "micro":100000  },
      "unittest"     : { "testfile":"polymers.test", "stable":False }
      }
    }

# generate json file
print >> open('polyelectrolytes.json', 'w+'), json.dumps(d, indent=4)

# generate molecule file with four atoms
print >> open('polymers.aam', 'w+'), """4
MM  0   0.00   0.00   0.00    1.0   1  3.0
MM  1   7.60   0.00   0.00    1.0   1  3.0
MM  2   0.00   7.60   0.00    1.0   1  3.0
MM  3   7.60   7.60   0.00    1.0   1  3.0"""

exe='./polyelectrolytes'
rc=1
if ( os.access( exe, os.X_OK )):
  rc = call( [exe] )
sys.exit( rc )

