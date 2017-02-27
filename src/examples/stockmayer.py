#!/usr/bin/env python
from __future__ import print_function
import json, sys, os
from subprocess import call, check_output
from shutil import copyfile

pfx = os.path.join( os.path.dirname(__file__), "stockmayer")
try:
  copyfile( pfx+'.state', 'state' )
  copyfile( pfx+'.test', 'stockmayer.test' )
except: pass

d = {
  "atomlist" : {
    "dip" : { "sigma":2.8863, "eps":1.9700, "mu":"0 0 1.652", "dp":0.7, "dprot":270, "alpha":"0 0 0 0 0 0" }
  },

  "moleculelist" : {
    "sol" : { "atoms":"dip", "Ninit":300, "atomic":True }
  },

  "energy" : {
    "nonbonded" : {
       "epsr":1.0, "cutoff":9.9188, "eps_rf":140.0
    }
  },
  
    "analysis" : {
      "xtcfile" :   { "file": "stockmayer.xtc", "nstep":20 },
      "energyfile": { "file": "energy.dat", "nstep":20 },
      "pqrfile" :   { "file": "stockmayer.pqr" },
      "statefile" : { "file": "state" },
      "multipoleanalysis" : { "nstep":20, "cutoff":9, "diel_RF":True },
      "kirkwoodfactor" : { "nstep":20, "pairs" :
            [
               { "name1":"dip", "name2":"dip", "dim":3, "dr":0.1, "file":"kwfactor_dipdip.dat"}
            ]
          },
      "atomrdf" : { "nstep":10, "pairs" :
            [
               { "name1":"dip", "name2":"dip", "dim":3, "dr":0.1, "file":"rdf_dipdip.dat"}
            ]
          }
      },

  "moves" : {
    "atomtranslate" : { "sol" : { "peratom":True } },
    "atomrotate" : { "sol" : { "peratom":True } }
  },

  "system" : {
    "temperature"  : 315.8,
    "geometry"     : { "length" : 19.8377 },
    "mcloop"       : { "macro":10, "micro":100 },
    "unittest"     : { "testfile":"stockmayer.test", "stable":False }
  }
}

# generate json file
with open('stockmayer.json', 'w+') as f:
    f.write(json.dumps(d, indent=4))

exe='./stockmayer'
rc=1
if ( os.access( exe, os.X_OK )):
  rc = call( [exe] )
sys.exit( rc )
