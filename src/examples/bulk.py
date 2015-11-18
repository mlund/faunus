#!/usr/bin/env python
import numpy as np
import json, sys, os
from subprocess import call, check_output
from shutil import copyfile

pfx = os.path.join( os.path.dirname(__file__), "bulk")
try:
  copyfile( pfx+'.state', 'state' )
  copyfile( pfx+'.test', 'bulk.test' )
except: pass

d = {
    "atomlist" : {
      "Na" : { "q": 1.0, "sigma":3.33, "eps":0.01158968, "dp":1.0 },
      "Cl" : { "q":-1.0, "sigma":4.40, "eps":0.4184,     "dp":1.0 }
      },

    "moleculelist" : {
      "NaCl" : { "atoms":"Na Cl", "Ninit":1152, "atomic":True }
      },

    "energy" : {
      "nonbonded" : { "coulomb" : { "epsr":1, "cutoff":14 } }
      },

    "analysis" : {
      "virial" : { "nstep":10 }
      },

    "moves" : {
      "isobaric" : { "dp":0.0, "pressure":11, "prob":0.1 },
      "atomtranslate" : { 
        "NaCl" : { "peratom":True }
        }
      },

    "system" : {
      "temperature"  : 1100,
      "cuboid"       : { "len" : 80 },
      "mcloop"       : { "macro":5, "micro":40 },
      "unittest"     : { "testfile":"bulk.test", "stable":False },
      "atomlist"     : "bulk.json",
      "moleculelist" : "bulk.json"
      }
    }

# generate json file
print >> open('bulk.json', 'w+'), json.dumps(d, indent=4)

exe='./bulk'
rc=1
if ( os.access( exe, os.X_OK )):
  rc = call( [exe] )
sys.exit( rc )
