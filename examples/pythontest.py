#!/usr/bin/env python

import json
import numpy as np
from pyfaunus import *

# Dictionary defining input
d = {}
d['temperature'] = 298.15
d['geometry'] = { 'length': 50 }
d['atomlist'] = [
        { 'Na': dict( r=2.0, eps=0.05, q=1.0, tfe=1.0 ) },
        { 'Cl': dict( r=1.2, eps=0.05, q=-1.0, tfe=1.0 ) }
        ]
d['moleculelist'] = [
        { 'salt': dict(atoms=['Na','Cl'], atomic=True) }
        ]
d['insertmolecules'] = [
        { 'salt': dict( N=1 ) }
        ]

# Create a simulation Space
#   this will:
#   - initialize `atoms` from dictionary
#   - add particles / molecules
spc = Space()
spc.from_dict(d)

# Loop over atom types
for i in atoms:
    print("atom name and diameter:", i.name, i.sigma)

# Use a pair-potential
print('\nSASA from pair-potential:\n')
d['energy'] = [
        { 'nonbonded': {
            'default' : [
                { 'sasa': { 'molarity': 1.0, 'radius': 1.4 } }
                ]
            } }
        ]
pot = FunctorPotential( d['energy'][0]['nonbonded'] )
r = np.linspace(0,10,10)
u = np.array( [pot.energy( spc.p[0], spc.p[1], [0,0,i] ) for i in r] )
print( np.array([r,u]) )

# Now let's test a Hamiltonian
print('\nSASA from PowerSASA (Hamiltonian):\n')
d['energy'] = [ { 'sasa': { 'molarity': 1.0, 'radius': 1.4 } } ]
H = Hamiltonian(spc, d)
spc.p[0].pos = [0,0,0] # fix 1st particle in origin
c = Change()  # change object telling that a full energy calculation
c.all = True; # should be performed when calling `energy()`
for i in r:   # loop over particle-particle distances
    spc.p[1].pos = [0,0,i]
    print ( H.energy(c) )

