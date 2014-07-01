#!/usr/bin/env python

import json, sys, os
from subprocess import call

def mkinput():
  # create input file
  def fmt(key,val,info): return "%-25s %-20s # %s\n" % (key,val,info)
  print >> open('./grand.input', 'w+'),\
      fmt("loop_macrosteps",   10, "number of macro loops"),\
      fmt("loop_microsteps", 1000, "number of micro loops"),\
      fmt("cuboid_len",   120, "cubic box length (angstrom)"),\
      fmt("temperature", 300, "(kelvin)"),\
      fmt("epsilon_r",    80, "dielectric constant"),\
      fmt("cardinaux_eps", 4.1, "cardinaux LJ epsilon (kT)"),\
      fmt("cardinaux_alpha", 90, "cardinaux alpha (integer)"),\
      fmt("atomlist", "grand.json", "atom property file"),\
      fmt("tion1", "M",  "ion type 1"),\
      fmt("tion2", "Na", "ion type 2"),\
      fmt("tion3", "Cl", "ion type 3"),\
      fmt("nion1", 10, "initial number of ion type 1"),\
      fmt("nion2", 20, "initial number of ion type 2"),\
      fmt("nion3", 70, "initial number of ion type 3"),\
      fmt("saltbath_runfraction", 1.0, "chance of running GC move (1=100%)"),\
      fmt("mv_particle_runfraction", 0.5, "chance of translating particles (1=100%)"),\

  # create atom property file (OBS: eps in kJ/mol and activity in mol/l)
  print >> open('./grand.json', 'w+'), json.dumps(
      {'atomlist': {
        'Na': {'q': 1.0, 'r':2.0,  'dp':50, 'activity':0.05, 'eps':0.001},
        'Cl': {'q':-1.0, 'r':2.0,  'dp':50, 'activity':0.05, 'eps':0.001},
        'M':  {'q': 5.0, 'r':10.0, 'dp':10, 'activity':0.0,  'eps':10}
        }}, indent=2)

mkinput()
call(["./boj-grand"])

