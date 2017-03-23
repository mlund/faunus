#!/usr/bin/env python

from __future__ import print_function
import json, sys, os
from subprocess import call, check_output
from shutil import copyfile

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--norun', action='store_true')
args = parser.parse_args()

name = 'water2'
pfx = os.path.join(os.path.dirname(__file__), name)
try:
    copyfile(pfx + '.state', 'state')
    copyfile(pfx + '.test', name + '.test')
except:
    pass

d = {
    'system': {
	"temperature"  : 300,
        'geometry': {'length': 18.67},
        'mcloop': {'macro': 10, 'micro': 200},
        'unittest': {'testfile': 'water2.test', 'stable': False}
    },
    'energy': {
        'nonbonded': {
            'coulombtype':"reactionfield" ,'epsr': 1.0, 'cutoff': 9.0, 'eps_rf':78.5, 'order':2, 'alpha':0.2,
            'ewald': {'eps_surf': 1e11, 'cutoff': 9, 'alpha': 0.2, "cutoffK":3, "spherical_sum":True, "update_frequency":216},
            'cutoff_g2g': 10
        }
    },
    'atomlist': {
        'OW': {'q': -0.8476, 'sigma': 3.165491652452016, 'eps': 0.650299305951160},
        'HW': {'q': 0.4238, 'sigma': 0.0, 'eps': 0.0}
    },
    'moleculelist': {
        'water': {'structure': 'water2.aam', 'Ninit': 216, 'insdir': '1 1 1'}
    },
    'moves': {
        'moltransrot': {
            'water': {'dp': 0.5, 'dprot': 0.5, 'dir': '1 1 1', 'permol': True, 'prob': 1.0}
        },
        'isobaric': {'dp': 0.1, 'pressure': 40.0906, 'prop': 1.0}
    },
    'analysis' : {
        'xtcfile' :   { 'file': 'water2.xtc', 'nstep':20 },
        'pqrfile' :   { 'file': 'water2.pqr' },
        'statefile' : { 'file': 'state' },
        "energyfile": { "file": "energy.dat", "nstep":20 },
        "multipoleanalysis" : { "nstep":20, "cutoff":9.0, "dielectric":"reactionfield", 'eps_rf':78.5 },
        'widommolecule' : dict(nstep=20, ninsert=10, molecule="water"),
        'sofq' :      dict(nstep=20, qmin=2, qmax=10, dq=0.5, mollist=["water"], file='debye.dat'),
        'atomrdf' : {
            'nstep':20, 'pairs' :
                [
                    { 'name1':'OW', 'name2':'OW', 'dim':3, 'file':'rdf_owow.dat', 'dr':0.05  }
                ]
            },
        'molrdf' : {
            'nstep':20, 'pairs' :
                [
                    { 'name1':'water', 'name2':'water', 'dim':3, 'file':'rdf_ww.dat', 'dr':0.05  }
                ]
            }
        }
}

# generate json file
with open(name + '.json', 'w+') as f:
    f.write(json.dumps(d, indent=4))

f = open('water2.aam', 'w')
f.write('''3
OW    1    2.30    6.28    1.13 -0.8476  15.99940 1.582745826226008
HW    2    1.37    6.26    1.50  0.4238   1.00800 1
HW    3    2.31    5.89    0.21  0.4238   1.00800 1
''')
f.close()

if args.norun == False:
    exe = './' + name
    rc = 1
    if (os.access(exe, os.X_OK)):
        rc = call([exe])
    sys.exit(rc)
