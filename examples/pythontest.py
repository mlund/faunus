import json
import unittest
import pyfaunus as pyf

# define atom types
d = { 'atomlist': [
    { 'Na+':  dict( r=2.0, eps=0.05, q=1.0 ) },
    { 'Cl-':  dict( r=1.2, q=-1.0 ) }
    ] }

pyf.atoms.from_dict(d)

a = pyf.atoms[0].p
b = pyf.atoms[1].p

b.charge = 10

# loop over atom types
for i in pyf.atoms:
    print("atom name and diameter:", i.name, i.sigma)

print("charge 2nd atomtype = ", pyf.atoms[1].p.charge)

#pot = pyf.FunctorPotential(d)
