#!/usr/bin/env python
from __future__ import print_function
import json, sys, os
from subprocess import call
from shutil import copyfile

rho=0.291    # reduced density
T=0.18       # reduced temperature
N=1000       # number of particles

L=(N/rho)**(0.5)
eps=1/T

j = {
        "atomlist" : { "CS" : { "r":0.5, "dp":0.9 } },
        "moleculelist" : {
            "mymol" : { "atoms":"CS", "atomic":True, "Ninit" : N, "insdir":"1 1 0" }
            },
        "energy" : {
            "nonbonded" : {
                "coreshell" : { "core_radius":1.0, "shell_radius":2.5, "epsilon" : eps }
                }
            },
        "moves" : {
            "atomtranslate" : {
                "mymol" : { "peratom":True, "dir":"1 1 0" }
                }
            },
        "system" : {
            "cuboid" : { "len" : L }
            },
        "analysis" : {
            "pqr" :       { "file": "stripes.pqr" },
            "statefile" : { "file": "state" }
            }
        }

with open("stripes.json", "w+") as f:
    f.write(json.dumps(j, indent=4))

exe="./stripes"
if ( os.access( exe, os.X_OK )):
    rc = call( [exe] )
  sys.exit( rc )
