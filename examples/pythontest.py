#!/usr/bin/env python

import json
import unittest
import numpy as np
from math import pi
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

# Analysis
analysis = Analysis(spc, H, d)


# Geometry

class TestGeometry(unittest.TestCase):

    def test_cuboid(self):
        geo = Chameleon( dict(type="cuboid", length=[2,3,4]) )
        V = geo.getVolume();
        self.assertAlmostEqual(V, 2*3*4, msg="volume")

        a = geo.boundary( [1.1, 1.5, -2.001] );
        self.assertAlmostEqual(a[0], -0.9)
        self.assertAlmostEqual(a[1], 1.5)
        self.assertAlmostEqual(a[2], 1.999)
        geo.setVolume(123.4);
        self.assertAlmostEqual( geo.getVolume(), 123.4);

    def test_sphere(self):
        r = 15.0
        geo = Chameleon( dict(type="sphere", radius=r) )
        V = geo.getVolume();
        self.assertAlmostEqual(V, 4*pi/3*r**3, msg="check volume")
        self.assertEqual( geo.collision([0,r+0.001,0]), True, msg="collision")
        self.assertEqual( geo.collision([0,0,r-0.001]), False, msg="no collision")
        geo.setVolume(123.4);
        self.assertAlmostEqual( geo.getVolume(), 123.4);

if __name__ == '__main__':
    unittest.main()

