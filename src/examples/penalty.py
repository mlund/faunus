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
      "xy-position": {
        "molecule":"myparticle", "f0":f0, "scale":0.7, "update":2e4, "bw1":0.1, "bw2":0.1, "lo1":-2, "hi1":2, "lo2":-2, "hi2":2 
        } 
    }
  },
  
  "moves" : { "atomtranslate" : { "myparticle" : { "peratom":True, "prob":1 } } },

  "system" : {
    "cuboid"   : { "len":4 },
    "mcloop"   : { "macro":10, "micro":micro },
    "unittest" : { "testfile":"penalty.test", "stable":False }
  }
}
  print >> open('penalty.json', 'w+'), json.dumps(j, indent=4)

exe  = './penalty'
if os.path.exists('./pf_penalty'):
  os.remove('./pf_penalty')
else:
  print("Sorry, I can not remove pf_penalty.")

if ( os.access( exe, os.X_OK )):

  f0=0
  micro=50000
  mkinput()
  rc = call( [exe] )
  f0=0.5
  micro=100000
  os.remove('./pf_penalty') 
  mkinput()
  rc = call( [exe] )
  f0=0
  micro=50000
  mkinput()
  rc = call( [exe] ) 

sys.exit( rc )

