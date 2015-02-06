#!/usr/bin/env python

import json, sys, os
from subprocess import call
from shutil import copyfile

rho=0.291    # reduced density
T=0.18       # reduced temperature
N=1000       # number of particles
L=(N/rho)**(0.5)
eps=1/T

js = {}
js['system'] = {}
js['system']['cuboid']       = { 'len':L }
js['system']['atomlist']     = 'stripes.json'
js['system']['moleculelist'] = 'stripes.json'

js['atomlist'] = {} 
js['atomlist']['CS'] = { 'r':0.5, 'dp':0.9 }

js['moleculelist'] = {} 
js['moleculelist']['mymol'] = { 'atoms':'CS', 'atomic':True, 'Ninit':N, 'insdir':'1 1 0' }

js['energy'] = {} 
js['energy']['nonbonded'] = { 'coreshell' : { 'core_radius':1.0, 'shell_radius':2.5, 'epsilon':eps }  }

js['moves'] = {}
js['moves']['atomtranslate'] = { 'mymol' : { 'peratom':True, 'dir':'1 1 0' } }

print >> open('stripes.json', 'w+'), json.dumps(js, indent=4)

exe='./stripes'
if ( os.access( exe, os.X_OK )):
  rc = call( [exe] )
  sys.exit( rc )
