#!/usr/bin/env python
from __future__ import print_function
import json, sys, os
from subprocess import call, check_output
from shutil import copyfile

pfx = os.path.join( os.path.dirname(__file__), "spheresurface")
try:
  copyfile( pfx+'.state', 'state' )
  copyfile( pfx+'.test', 'spheresurface.test' )
except: pass

d = {
    "atomlist" : {
      "Pos" : { "q": 1.0, "sigma":1.0, "eps":0.0, "dp":0.2 },
      "Neg" : { "q":-1.0, "sigma":1.0, "eps":0.0, "dp":0.2 }
      },

    "moleculelist" : {
      "Charges" : { "atoms":"Pos Neg", "Ninit":100, "atomic":True }
      },

    "energy" : {
        "nonbonded" : { "coulomb" : { "epsr":72.0 } }
      },

    "analysis" : {
      "xtcfile" :   { "file": "spheresurface.xtc", "nstep":20 },
      "energyfile": { "file": "energy.dat", "nstep":20 },
      "pqrfile" :   { "file": "spheresurface.pqr" },
      "statefile" : { "file": "state" },
      "atomrdf" : { "nstep":10, "pairs" :
            [
               { "name1":"Pos", "name2":"Neg", "dim":2, "dr":0.1, "file":"rdf_PosNeg.dat", "Rhyper":10.0},
               { "name1":"Pos", "name2":"Pos", "dim":2, "dr":0.1, "file":"rdf_PosPos.dat", "Rhyper":10.0}
            ]
          }
      },

    "moves" : {
      "atomictranslation2D" : { 
          "Charges" : { "peratom":True, "prob": 1.0 }
        }
      },

    "system" : {
      "temperature"  : 300,
      "geometry"     : { "radius" : 10.0 },
      "mcloop"       : { "macro":10, "micro":100 },
      "unittest"     : { "testfile":"spheresurface.test", "stable":False },
      "atomlist"     : "spheresurface.json",
      "moleculelist" : "spheresurface.json"
      }
    }

# generate json file
with open('spheresurface.json', 'w+') as f:
    f.write(json.dumps(d, indent=4))

exe='./spheresurface'
rc=1
if ( os.access( exe, os.X_OK )):
  rc = call( [exe] )
sys.exit( rc )
