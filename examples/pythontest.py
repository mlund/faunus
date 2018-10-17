import json
import pyfaunus as pyf

inputstr = '''
{
  "atomlist" : [
    { "A": { "r":1.1 } },
    { "B": { "activity":0.2, "eps":0.05, "dp":9.8, "dprot":3.14, "mw":1.1, "tfe":0.98, "tension":0.023 } }
  ]
}
'''
# define atom types
d = {'atomlist': [
  dict(A = dict( r=2.0, eps=0.05 )),
  dict(B = dict( r=1.2 )) ]}

pyf.atoms.from_dict(d)

# loop over all atom types
for i in pyf.atoms:
    print("atom name and diameter:", i.name, i.sigma )

pot = pyf.FunctorPotential(d)
