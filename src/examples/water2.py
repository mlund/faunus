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
        'geometry': {'length': 32.00},
        'mcloop': {'macro': 10, 'micro': 50},
        'unittest': {'testfile': 'water2.test', 'stable': False}
    },
    'energy': {
        'nonbonded': {
            'coulombtype':"fanourgakis" ,'epsr': 1.0, 'cutoff': 14.8, 'eps_rf':78.5, 'order':2, 'alpha':0.2,
            'ewald': {'eps_surf': 1e11, 'cutoff': 14.8, 'alpha': 0.21, "cutoffK":7.3, "spherical_sum":True, "update_frequency":1000},
            'cutoff_g2g': 15.0
        }
    },
    'atomlist': {
        'OW': {'q': -0.8476, 'sigma': 3.165491652452016, 'eps': 0.650299305951160},
        'HW': {'q': 0.4238, 'sigma': 0.0, 'eps': 0.0}
    },
    'moleculelist': {
        'water': {'structure': 'water2.aam', 'Ninit': 1000, 'insdir': '1 1 1'}
    },
    'moves': {
        'moltransrot': {
            'water': {'dp': 0.4, 'dprot': 0.4, 'dir': '1 1 1', 'permol': True, 'prob': 1.0}
        },
        'isobaric': {'dp': 0.03, 'pressure': 40.0906, 'prop': 1.0}
    },
    'analysis' : {
        'chargemultipole' : {'nstep':20, 'mollist':['water'] },
        'xtcfile' :   { 'file': 'water2.xtc', 'nstep':20 },
        'pqrfile' :   { 'file': 'water2.pqr' },
        'statefile' : { 'file': 'state' },
        "energyfile": { "file": "energy.dat", "nstep":20 },
        "multipoleanalysis" : { "nstep":20, "cutoff":14.8, "dielectric":"tinfoil", 'eps_rf':78.5, 'pairs' :
                [
                    { 'name1':'water', 'name2':'water', 'dim':3, 'file':'rdf_ww.dat', 'file2':'mucorr_ww.dat', 'dr':0.05  }
                ] 
		},
        'widommolecule' : dict(nstep=20, ninsert=10, molecule="water"),
        'scatter' :      dict(nstep=20, qmin=2, qmax=10, dq=0.5, mollist=["water"], file='debye.dat'),
        'atomrdf' : {
            'nstep':20, 'pairs' :
                [
                    { 'name1':'OW', 'name2':'OW', 'dim':3, 'file':'rdf_owow.dat', 'dr':0.05  }
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
