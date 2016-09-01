#!/usr/bin/env python

#SBATCH -N 1
#SBATCH --tasks-per-node=2
#SBATCH -t 00:10:00
#SBATCH -J test
#SBATCH --qos=test

from __future__ import print_function
import numpy as np
import json, sys, os
from subprocess import call
from shutil import copyfile

name='manybody'
exe='./'+name
pfx = os.path.join( os.path.dirname(__file__), name)
try:
  copyfile( pfx+'.state', 'state' )
  copyfile( pfx+'.test', name+'.test' )
except: pass

def mkinput():
  d = {
    "energy" : {
      "hydrophobicsasa" : {},
      "nonbonded" : {
        "cutoff_g2g" : cut_g2g,
        "ljsimple" : { "eps":0.05, "cutoff":cut_i2i },
        "coulomb" : { "epsr" : 78.7, "ionicstrength":salt, "cutoff":cut_i2i }
      }
    },

    "atomlist" : {
      "H3PO4":  { "q":0,  "r":2.0 },
      "H2PO4":  { "q":-1, "r":2.0 },
      "HPO4" :  { "q":-2, "r":2.0 },
      "PO4"  :  { "q":-3, "r":2.0 },
      "BPTI" :  { "q":7.3, "r":12.29 },
      "Na"   :  { "q": 1, "r":1.9, "mw":22.99 },
      "Cl"   :  { "q":-1, "r":1.7, "mw":35.45 },
      "I"    :  { "q":-1, "r":2.0, "mw":1 },
      "SCN"  :  { "q":-1, "r":2.0, "mw":1 },
      "ASP"  :  { "q":-1, "r":3.6, "mw":110 },
      "HASP" :  { "q":0,  "r":3.6, "mw":110 },
      "LASP" :  { "q":2,  "r":3.6, "mw":110 },
      "CTR"  :  { "q":-1, "r":2.0, "mw":16 },
      "HCTR" :  { "q":0,  "r":2.0, "mw":16 },
      "GLU"  :  { "q":-1, "r":3.8, "mw":122 },
      "HGLU" :  { "q":0,  "r":3.8, "mw":122 },
      "LGLU" :  { "q":2,  "r":3.8, "mw":122 },
      "HIS"  :  { "q":0,  "r":3.9, "mw":130 },
      "HHIS" :  { "q":1,  "r":3.9, "mw":130 },
      "NTR"  :  { "q":0,  "r":2.0, "mw":14 },
      "HNTR" :  { "q":1,  "r":2.0, "mw":14 },
      "TYR"  :  { "q":-1, "r":4.1, "mw":154 },
      "HTYR" :  { "q":0,  "r":4.1, "mw":154 },
      "LYS"  :  { "q":0,  "r":3.7, "mw":116 },
      "HLYS" :  { "q":1,  "r":3.7, "mw":116 },
      "CYS"  :  { "q":-1, "r":3.6, "mw":103 },
      "HCYS" :  { "q":0,  "r":3.6, "mw":103 },
      "CYb"  :  { "q":0,  "r":3.6, "mw":103 },
      "ARG"  :  { "q":0,  "r":4.0, "mw":144 },
      "HARG" :  { "q":1,  "r":4.0, "mw":144 },
      "ALA"  :  { "q":0,  "r":3.1, "mw":66 },
      "ILE"  :  { "q":0,  "r":3.6, "mw":102 },
      "LEU"  :  { "q":0,  "r":3.6, "mw":102 },
      "MET"  :  { "q":0,  "r":3.8, "mw":122 },
      "PHE"  :  { "q":0,  "r":3.9, "mw":138 },
      "PRO"  :  { "q":0,  "r":3.4, "mw":90 },
      "TRP"  :  { "q":0,  "r":4.3, "mw":176 },
      "VAL"  :  { "q":0,  "r":3.4, "mw":90 },
      "SER"  :  { "q":0,  "r":3.3, "mw":82 },
      "THR"  :  { "q":0,  "r":3.5, "mw":94 },
      "ASN"  :  { "q":0,  "r":3.6, "mw":108 },
      "GLN"  :  { "q":0,  "r":3.8, "mw":120 },
      "GLY"  :  { "q":0,  "r":2.9, "mw":54 }
    },

    "processes" : {
      "H-Asp" : { "bound":"HASP" , "free":"ASP" , "pKd":4.0  , "pX":pH },
      "H-Ctr" : { "bound":"HCTR" , "free":"CTR" , "pKd":2.6  , "pX":pH },
      "H-Glu" : { "bound":"HGLU" , "free":"GLU" , "pKd":4.4  , "pX":pH },
      "H-His" : { "bound":"HHIS" , "free":"HIS" , "pKd":6.3  , "pX":pH },
      "H-Arg" : { "bound":"HARG" , "free":"ARG" , "pKd":12.0 , "pX":pH },
      "H-Ntr" : { "bound":"HNTR" , "free":"NTR" , "pKd":7.5  , "pX":pH },
      "H-Cys" : { "bound":"HCYS" , "free":"CYS" , "pKd":10.8 , "pX":pH },
      "H-Tyr" : { "bound":"HTYR" , "free":"TYR" , "pKd":9.6  , "pX":pH },
      "H-Lys" : { "bound":"HLYS" , "free":"LYS" , "pKd":10.4 , "pX":pH },
      "K1"    : { "bound":"H3PO4", "free":"H2PO4","pKd":2.12,  "pX":pH },
      "K2"    : { "bound":"H2PO4", "free":"HPO4", "pKd":7.21,  "pX":pH },
      "K3"    : { "bound":"HPO4",  "free":"PO4",  "pKd":12.67, "pX":pH }
    },

    "moleculelist": {
      "protein":  { "structure":"manybody.bpti", "Ninit":Np }
    },

    "moves" : {
       "titrate" : { "prob":0.2, "processfile":"manybody.json" },
       "moltransrot" : {
         "protein" : { "dp":dp, "dprot":dprot, "prob":1.0, "permol":True } 
       } 
    },

    "analysis" : {
      "scatter" : { "qmin":0.05, "qmax":0.4, "dq":0.005, "cutoff":1e10 }
    },

    "system" : {
      "temperature" : 298.15,
      "cuboid"   : { "len" : boxlen },
      "mcloop"   : { "macro":10, "micro":micro },
      "unittest" : { "testfile":"manybody.test", "stable":False }
    }
  }
  with open(name+'.json', 'w+') as f:
      f.write(json.dumps(d, indent=4))

Zp    = 7.0   # Protein charge (approximate, for counter ion calc.)
Rp    = 20.   # Protein max radius used for cutoff determination
Np    = 20    # Number of protein molecules in simulation
Cs    = 0.020 # 1:1 salt concentration (mol/l)
dp    = 60    # translational displacement (angstrom)
dprot = 6     # rotational displacement

# Flow control
eqrun   = True
prodrun = False
copy    = False

# Protein concentration [mol/l]
for Cp in [0.004]:

  V       = Np / (1e-27*6.022e23*Cp)   # simulation volume
  boxlen  = V**(1/3.)                  # box length (angstrom)
  salt    = 0.5*Zp*Cp + Cs             # counter ion contribution to ionic strength
  D       = 3.04/salt**0.5             # debye length
  cut_i2i = 3*D                        # residue-residue cutoff
  cut_g2g = 2*Rp+cut_i2i               # protein-protein COM cutoff

  for pH in [4.1]:
    prefix='pH'+str(pH)+'-Cp'+str(Cp)
    print(prefix)

    if eqrun==True:
      print('Equilibration run...(state file deleted)')
      try:
        os.remove('state')
      except: pass
      micro=100
      mkinput()
      call( ['mpiexec', exe] )

    if prodrun==True:
      print('Production run...')
      micro=1000
      mkinput()
      call( ['mpiexec', exe] )

    if copy==True:
      os.rename('state', prefix+'.state')
      os.rename('debye.dat', prefix+'.debye')
      os.rename('rdf_p2p.dat', prefix+'.rdf')
      os.rename('confout.pqr', prefix+'.pqr')
      os.rename('avgcharge.pqr', prefix+'.avgcharge.pqr')
      os.rename('manybody.input', prefix+'.input')

