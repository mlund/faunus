import json, sys, os
from subprocess import call
from shutil import copyfile

def mkinput():
  # create input file
  def fmt(key,val,info): return "%-25s %-20s # %s\n" % (key,val,info)
  print >> open('./grand.input', 'w+'),\
      fmt("temperature", 300, "(kelvin)"),\
      fmt("epsilon_r", 80, "dielectric constant"),\
      fmt("sphere_radius", 80, "radius of container (angstrom)"),\
      fmt("atomlist", "grand.json", "atom property file"),\
      fmt("loop_macrosteps", 10, "number of macro loops"),\
      fmt("loop_microsteps", 20000, "number of micro loops"),\
      fmt("tion1", "Na", "ion type 1"),\
      fmt("tion2", "Cl", "ion type 2"),\
      fmt("nion1", 20, "initial number of ion type 1"),\
      fmt("nion2", 20, "initial number of ion type 2"),\
      fmt("saltbath_runfraction", 1.0, "chance of running GC move (1=100%)"),\
      fmt("mv_particle_runfraction", 0.01, "chance of translating particles (1=100%)"),\
      fmt("test_stable", "no", "create (yes) or check test file (no)"),\
      fmt("test_file", "grand.test", "name of test file to create or load")

  # create json file
  print >> open('./grand.json', 'w+'), json.dumps(
      {'atomlist': {
        'Na': {'q': 1.0, 'r':2.0, 'dp':50, 'activity':0.05},
        'Cl': {'q':-1.0, 'r':2.0, 'dp':50, 'activity':0.05}
        }}, indent=2)

exe="./grand"
if (os.access(exe,os.X_OK)):
  copyfile("grand.state", "state")
  mkinput()
  rc = call([exe])
  sys.exit(rc)

