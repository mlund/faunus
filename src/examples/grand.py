import json, sys, os
from subprocess import call
from shutil import copyfile

def mkinput():
  js = {}

  # direct addition
  js['atomlist'] = {} 
  js['atomlist']['Na'] = { 'q': 1.0, 'r':2.0, 'dp':50, 'activity':0.05 }
  js['atomlist']['Cl'] = { 'q':-1.0, 'r':2.0, 'dp':50, 'activity':0.05 }

  js['moleculelist'] = {} 
  js['moleculelist']['salt'] = { 'atoms':'Na Cl', 'atomic':True, 'Ninit':20 }

  js['energy'] = {} 
  js['energy']['nonbonded'] = { 'coulomb' : { 'epsr':80 }  }

  js['moves'] = {}
  js['moves']['atomgc']        = { 'molecule':'salt' }
  js['moves']['atomtranslate'] = { 'salt' : { 'prob':0.01, 'permol':True } }

  # adding via reference
  sec = js['system'] = {}
  sec['temperature']  = 298.15
  sec['sphere']       = { 'radius':80 }
  sec['mcloop']       = { 'macro':10, 'micro':100000 }
  sec['coulomb']      = { 'epsr':80 }
  sec['unittest']     = { 'testfile':'grand.test', 'stable':False }
  sec['atomlist']     = 'grand.json'
  sec['moleculelist'] = 'grand.json'

  print >> open('grand.json', 'w+'), json.dumps(js, indent=4)

exe='./grand'
if ( os.access( exe, os.X_OK )):
  copyfile( 'grand.state', 'state' )
  mkinput()
  rc = call( [exe] )
  sys.exit( rc )

