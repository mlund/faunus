#!/usr/bin/env python
import json, sys, os
from subprocess import call, check_output
from shutil import copyfile

pfx = os.path.join( os.path.dirname(__file__), "slitpolymer")
try:
  copyfile( pfx+'.state', 'state' )
  copyfile( pfx+'.test', 'slitpolymer.test' )
except: pass

d = {
    "atomlist" : {
      "MM" : { "q": 1.0, "r":3.0, "dp":10 }
      },

    "moleculelist" : {
      "polymer" : { "structure":"slitpolymer.aam", "Ninit":1 }
      },

    "analysis" : {
      "polymershape" : { "nstep":10, "mollist":[ "polymer" ] }
      },

    "energy" : {
      "gouychapman" : { "phi0":-2 },
      "nonbonded" : { "coulomb" : { "epsr" : 78.9, "debyelength":11.6492873993050 } }
      },

    "moves" : {
      "moltransrot" : {
        "polymer" : { "dp":100, "dprot":6, "permol":True } 
        },
      "crankshaft" : {
        "polymer" : { "dp":6, "maxlen":6, "permol":True }
        },
      "pivot" : {
        "polymer" : { "dp":6, "maxlen":6, "permol":True }
        },
      "reptate" : {
        "polymer" : { "bondlength":4.9, "permol":True }
        }
      },

    "system" : {
      "temperature"  : 298,
      "cuboid"       : { "xyzlen" : "300 300 154" },
      "mcloop"       : { "macro":10, "micro":1000000 },
      "unittest"     : { "testfile":"slitpolymer.test", "stable":False }
      }
    }

# generate json file
print >> open('slitpolymer.json', 'w+'), json.dumps(d, indent=4)

f = open('slitpolymer.aam', 'w')
f.write('10\n\
    MM 1 -29.023 -32.227 0 0.82686 1 2.45\n\
    MM 2 -24.123 -32.227 0 -0.99969 1 2.45\n\
    MM 3 -19.223 -32.227 0 0 1 2.45\n\
    MM 4 -14.323 -32.227 0 0.098001 1 2.45\n\
    MM 5 -9.4231 -32.227 0 0 1 2.45\n\
    MM 6 -4.5231 -32.227 0 0.99884 1 2.45\n\
    MM 7 0.37692 -32.227 0 0.99998 1 2.45\n\
    MM 8 5.2769 -32.227 0 0.040585 1 2.45\n\
    MM 9 5.2769 -27.327 0 0.063099 1 2.45\n\
    MM 10 5.2769 -22.427 0 0 1 2.45\n')
f.close()

exe='./slitpolymer'
rc=1
if (os.access( exe, os.X_OK )):
  rc = call([exe])
sys.exit(rc)

