#!/usr/bin/env python

from __future__ import print_function
import json, sys, os
from subprocess import call, check_output
from shutil import copyfile

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--norun', action='store_true')
args = parser.parse_args()

name = 'gcmol'
pfx = os.path.join(os.path.dirname(__file__), name)
try:
    copyfile(pfx + '.state', 'state')
    copyfile(pfx + '.test', name + '.test')
except:
    pass

d = {
    'system': {
        'cuboid': {'len': 19.0288},
        'mcloop': {'macro': 10, 'micro': 100000},
        'unittest': {'testfile': 'gcmol.test', 'stable': False}
    },

    'energy': {
        'nonbonded': {
            'coulomb': {'epsr': 1.0, 'cutoff': 9.0}
        }
    },

    'atomlist': {
        'OW': {'q': -0.8476, 'sigma': 3.2, 'eps': 0.65},
        'HW': {'q': 0.4238, 'sigma': 0.0, 'eps': 0}
    },

    'moleculelist': {
        'water': {'activity': 0.0009386, 'structure': 'water.aam', 'Ninit': 200}
    },

    'analysis': {
        'pqrfile' :   { 'file': 'confout.pqr' },
        'statefile' : { 'file': 'state' }
        },

    'moves': {
        'gc': {
            'moleculecombinations': {
                'water': {'molecules': 'water'}
            },
        },
        'moltransrot': {
            'water': {'dp': 0.5, 'dprot': 0.5, 'prob': 0.001, 'permol': True}
        }
    }
}

# generate json file
with open(name + '.json', 'w+') as f:
    f.write(json.dumps(d, indent=4))

# Make SPC water molecule
f = open('water.aam', 'w')
f.write('''3
OW    1    2.30    6.28    1.13 -0.8476  15.99940 1.6
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
