#!/usr/bin/env python
from __future__ import print_function
import json, sys, os
from subprocess import call, check_output
from shutil import copyfile

pfx = os.path.join( os.path.dirname(__file__), "capparticles")
try:
  copyfile( pfx+'.state', 'state' )
  copyfile( pfx+'.test', 'capparticles.test' )
except: pass

d = {
    "atomlist" : {
      "Cap" : { "q": 0.0, "sigma":1.00, "dp":1.0, "dprot":360.0, "cap_center":0.5, "charge_position":0.0, "cap_radius":0.5 },
      "Sph" : { "q": 0.0, "sigma":1.00, "dp":1.0, "dprot":0.0, "cap_center":0.0, "charge_position":0.0, "cap_radius":0.0 }
      },

    "moleculelist" : {
      "Capsphere" : { "atoms":"Cap", "Ninit":1, "atomic":True },
      "Sphere" : { "atoms":"Sph", "Ninit":10, "atomic":True }
      },

    "energy" : {
        "nonbonded" : {'coulombtype':"none" ,'epsr': 1.0, 'cutoff': 10, 'order':2, 'alpha':0.0 }
      },

    "analysis" : {
      "povray" :   { "file": "some.pov" },
      "xtcfile" :   { "file": "capparticles.xtc", "nstep":100000000 },
      "energyfile": { "file": "energy.dat", "nstep":20 },
      "pqrfile" :   { "file": "capparticles.pqr" },
      "statefile" : { "file": "state" },
      "atomrdf" : { "nstep":10, "pairs" :
            [
               { "name1":"Cap", "name2":"Cap", "dim":3, "dr":0.01, "file":"rdf_capcap.dat"},
               { "name1":"Cap", "name2":"Sph", "dim":3, "dr":0.01, "file":"rdf_capsph.dat"},
               { "name1":"Sph", "name2":"Sph", "dim":3, "dr":0.01, "file":"rdf_sphsph.dat"}
            ]
          },
      "capanalysis" : { "nstep":10, "pairs" :
            [
               { "name1":"Cap", "name2":"Cap", "dim":3, "dr":0.05, "file":"corr_capcap.dat"}
            ]
          }
      },

  "moves" : {
    "atomtranslate" : { 
			"Capsphere" : { "peratom":True, "prob":1.0 }, 
		        "Sphere" : { "peratom":True, "prob":1.0 } 
		        },
    "atomrotate" : { "Capsphere" : { "peratom":True, "prob":1.0 }  }
  },

    "system" : {
      "temperature"  : 300,
      "geometry"     : { "length" : 3.0 },
      "mcloop"       : { "macro":10, "micro":10000 },
      "unittest"     : { "testfile":"capparticles.test", "stable":False },
      "atomlist"     : "capparticles.json",
      "moleculelist" : "capparticles.json"
      }
    }

# generate json file
with open('capparticles.json', 'w+') as f:
    f.write(json.dumps(d, indent=4))

exe='./capparticles'
rc=1
if ( os.access( exe, os.X_OK )):
  rc = call( [exe] )
sys.exit( rc )
