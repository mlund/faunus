#!/usr/bin/env python
import json, sys, os
from subprocess import call
from shutil import copyfile

def mkinput():
  global pH, micro
  j = {
      "processes": {
        "M-910": { "pKd": 5,   "pX": pH, "bound": "M-09", "free": "M-10" },
        "M-89":  { "pKd": 4.5, "pX": pH, "bound": "M-08", "free": "M-09" },
        "M-78":  { "pKd": 4,   "pX": pH, "bound": "M-07", "free": "M-08" },
        "M-67":  { "pKd": 3.5, "pX": pH, "bound": "M-06", "free": "M-07" },
        "M-56":  { "pKd": 3,   "pX": pH, "bound": "M-05", "free": "M-06" },
        "M-45":  { "pKd": 2.5, "pX": pH, "bound": "M-04", "free": "M-05" },
        "M-34":  { "pKd": 2,   "pX": pH, "bound": "M-03", "free": "M-04" },
        "M-23":  { "pKd": 1.5, "pX": pH, "bound": "M-02", "free": "M-03" },
        "M-12":  { "pKd": 1,   "pX": pH, "bound": "M-01", "free": "M-02" },
        "M-01":  { "pKd": 0.5, "pX": pH, "bound": "M+-0", "free": "M-01" },
        "M+10":  { "pKd": -0.5,"pX": pH, "bound": "M+01", "free": "M+-0" },
        "M+21":  { "pKd": -1,  "pX": pH, "bound": "M+02", "free": "M+01" },
        "M+32":  { "pKd": -1.5,"pX": pH, "bound": "M+03", "free": "M+02" },
        "M+43":  { "pKd": -2,  "pX": pH, "bound": "M+04", "free": "M+03" },
        "M+54":  { "pKd": -2.5,"pX": pH, "bound": "M+05", "free": "M+04" },
        "M+65":  { "pKd": -3,  "pX": pH, "bound": "M+06", "free": "M+05" },
        "M+76":  { "pKd": -3.5,"pX": pH, "bound": "M+07", "free": "M+06" },
        "M+87":  { "pKd": -4,  "pX": pH, "bound": "M+08", "free": "M+07" },
        "M+98":  { "pKd": -4.5,"pX": pH, "bound": "M+09", "free": "M+08" },
        "M+109": { "pKd": -5,  "pX": pH, "bound": "M+10", "free": "M+09" }
        },

      "energy": {
        "eqstate": { "processfile": "gctit.json" },
        "nonbonded": {
          "coulomb": { "epsr": 80 },
          "ljsimple": { "eps": 0.05 }, "coulomb": { "epsr": 80 }
          }
        },

      "system": {
        "temperature": 298.15,
        "cuboid": { "len": 202.5 },
        "unittest": { "testfile": "gctit.test", "stable": False },
        "mcloop": { "macro": 10, "micro": micro },
        "sphere": { "radius": 100 }
        },

      "moves": {
        "gctit"         : { "molecule": "salt", "prob": 0.01 },
        "atomtranslate" : { "salt": { "prob": 0.1, "dp": 100 } },
        "moltransrot"   : { "protein":
          { "permol": True, "dp": 60, "prob": 1, "dir": "0 0 1", "dprot": 0 } }
        },

      "moleculelist": {
        "protein": { "Ninit": 2, "structure": "gctit_mol.aam", "insdir": "0 0 1" },
        "salt": { "Ninit": 50, "atomic": True, "atoms": "La Cl Cl Cl" }
        },

      "analysis" : {
        "widom" : { "nstep":20, "ninsert":15, "particles":["La", "Cl", "Cl", "Cl"] },
        "widomscaled" : { "nstep":20, "ninsert":15, "lB":7.0 },
        "molrdf" : {
            "nstep":20, "pairs" :
                [
                    dict(name1="protein", name2="protein", dim=3, file="rdf.dat", dr=0.1)
                ]
            }
        },

      "atomlist": {
          "La":   { "q": 3, "r": 2, "dp": 20, "activity": 0.001601 },
          "Cl":   { "q": -1,"r": 2, "dp": 50, "activity": 0.02276 },
          "Na":   { "q": 1, "r": 2, "dp": 50 },
          "ASP":  { "q": -1, "r": 2 },
          "HASP": { "q": 0, "r": 2 },
          "GLU":  { "q": -1, "r": 2 },
          "HGLU": { "q": 0, "r": 2 },
          "M+08": { "q": 8, "r": 5 },
          "M+09": { "q": 9, "r": 5 },
          "M+01": { "q": 1, "r": 5 },
          "M+02": { "q": 2, "r": 5 },
          "M+03": { "q": 3, "r": 5 },
          "M+04": { "q": 4, "r": 5 },
          "M+05": { "q": 5, "r": 5 },
          "M+06": { "q": 6, "r": 5 },
          "M+07": { "q": 7, "r": 5 },
          "M-08": { "q": -8, "r": 5 },
          "M-06": { "q": -6, "r": 5 },
          "M-07": { "q": -7, "r": 5 },
          "M-04": { "q": -4, "r": 5 },
          "M-05": { "q": -5, "r": 5 },
          "M-02": { "q": -2, "r": 5 },
          "M-03": { "q": -3, "r": 5 },
          "M-01": { "q": -1, "r": 5 },
          "M+10": { "q": 10, "r": 5 },
          "M-09": { "q": -9, "r": 5 },
          "M+-0": { "q": 0, "r": 5 },
          "M-10": { "q": -10, "r": 5 }
          }
      }
  with open('gctit.json', 'w+') as f:
    f.write(json.dumps(j, indent=4))

exe            = './gctit'
runeq          = True
runprod        = False

if ( os.access( exe, os.X_OK )):
  for activity in [ 0.011 ]: # mol/l
    for pH in [ 0.0 ]:

      if (runeq==True):
        try:
          os.remove('state')
        except OSError:
          pass
        micro = 10000
        mkinput()
        rc = call( [exe] ) 

      if (runprod==True):
        micro = 100000
        mkinput()
        rc = call( [exe] ) 

sys.exit( rc )

