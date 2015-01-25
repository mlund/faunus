import json, sys, os
from subprocess import call
from shutil import copyfile

def mkinput():

  # create json file
  print >> open('grand.json', 'w+'), json.dumps(
      {
        'atomlist': {
          'Na': {'q': 1.0, 'r':2.0, 'dp':50, 'activity':0.05},
          'Cl': {'q':-1.0, 'r':2.0, 'dp':50, 'activity':0.05}
          },
        'moleculelist': {
          'salt': {'atoms':'Na Cl', 'atomic':True, 'Ninit':20 }
          },
        'moves': {
          'atomtranslate' : {
            'salt' : { 'prop':0.01, 'permol':True }
            },
          'atomgc' : { 'molecule':'salt' }
          },
        'energy' : {
          'nonbonded' : { 'coulomb' : { 'epsr' : 80. } }
          },
        'system': {
          'temperature'  : 300,
          'sphere'       : { 'radius' : 80. },
          'coulomb'      : { 'epsr' : 80. },
          'mcloop'       : { 'macro':10, 'micro':100000 },
          'unittest'     : { 'testfile':'grand.test', 'stable':False },
          'atomlist'     : 'grand.json',
          'moleculelist' : 'grand.json'
          }
        }, 
      indent=4 )

exe="./grand"
if (os.access(exe,os.X_OK)):
  copyfile("grand.state", "state")
  mkinput()
  rc = call([exe])
  sys.exit(rc)

