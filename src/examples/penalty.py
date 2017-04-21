#!/usr/bin/env python

from __future__ import print_function
import json, sys, os
from subprocess import call
from shutil import copyfile

def mkinput():
  global f, micro
  j = {
  "atomlist" : { "A" : { "r":0, "dp":0.5 } },

  "moleculelist" : { "myparticle" : { "atoms":"A", "atomic":True, "Ninit":1 } },

  "energy" : { 
    "penalty" : { 
      "xyz": {
          "first":"myparticle", "f0":f0, "scale":0.5, "update":1e4, "bw1":0.1, "bw2":0.1, "lo1":-2, "hi1":2, "lo2":-2, "hi2":2, "dir":"1 1 0" 
        } 
    }
  },
  
  "moves" : { "atomtranslate" : { "myparticle" : { "peratom":True, "prob":1 } } },

  "system" : {
    "geometry" : { "length":4 },
    "mcloop"   : { "macro":10, "micro":micro },
    "unittest" : { "testfile":"penalty.test", "stable":False }
  }
}
  with open('penalty.json', 'w+') as f:
      f.write(json.dumps(j, indent=4))

exe  = './penalty'
if os.path.exists('./pf_penalty'):
  os.remove('./pf_penalty')

if ( os.access( exe, os.X_OK )):

  f0=0
  micro=50000
  mkinput()
  rc = call( [exe] )
  f0=0.1
  micro=100000
  mkinput()
  rc = call( [exe] )
  f0=0
  micro=50000
  mkinput()
  rc = call( [exe] ) 

sys.exit( rc )

