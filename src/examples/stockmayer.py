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
    "sol" : { "atoms":"dip", "Ninit":1000, "atomic":True }
  },

  "energy" : {
    "nonbonded" : {
       "coulombtype":"wolf" ,"epsr":1.0, "cutoff":14.4315, "eps_rf":140.0, "alpha":0.21,
       'ewald': { "cutoff":14.4315,'eps_surf':1e11, "alpha":0.14, "cutoffK":4, "spherical_sum":True, "update_frequency":1000, "isotropic_pbc":False}
    }
  },
  
    "analysis" : {
      "_xtcfile" :   { "file": "stockmayer.xtc", "nstep":20 },
      "energyfile": { "file": "energy.dat", "nstep":20 },
      "pqrfile" :   { "file": "stockmayer.pqr" },
      "statefile" : { "file": "state" },
      "multipoleanalysis" : { "nstep":20, "cutoff":14.4315, "dielectric":"reactionfield", "eps_rf":140.0, "kappa":0.21, 'pairs' :
            [
               { 'name1':'sol', 'name2':'sol', 'dim':3, 'file':'rdf_ss.dat', 'dr':0.1  }
            ] 
	  },
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
    "atomtranslate" : { "sol" : { "peratom":True, "prob":1.0 } },
    "atomrotate" : { "sol" : { "peratom":True, "prob":1.0 } }
  },

  "system" : {
    "temperature"  : 315.8,
    "geometry"     : { "length" : 29.6336 },
    "mcloop"       : { "macro":10, "micro":50 },
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
