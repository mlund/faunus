#!/usr/bin/env python

import json
import numpy as np
from pyfaunus import *

# Dictionary defining our input
d = {}
d['geometry'] = { 'length': 50 }
d['atomlist'] = [
        { 'Na': dict( r=2.0, eps=0.05, q=1.0, tfe=1.0 ) },
        { 'Cl': dict( r=1.2, eps=0.05, q=-1.0, tfe=0.5 ) }
        ]
d['moleculelist'] = [
        { 'salt': dict(atoms=['Na','Cl'], atomic=True) }
        ]
d['insertmolecules'] = [
        { 'salt': dict( N=1 ) }
        ]
d['energy'] = {
        'nonbonded': [
            { 'default' : [
                { 'sasa': { 'molarity': 1.0 } }
                ] }
            ] }

# Create a simulation Space
spc = Space()
spc.from_dict(d)

# Loop over atom types
for i in atoms:
    print("atom name and diameter:", i.name, i.sigma)

# Use a pair-potential
pot = FunctorPotential( d['energy']['nonbonded'][0] )
r = np.linspace(1,5,10)
u = np.array( [pot.energy( spc.p[0], spc.p[1], [0,0,i] ) for i in r] )
print( np.array([r,u]) )

