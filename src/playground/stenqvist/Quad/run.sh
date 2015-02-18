#!/bin/bash

echo '{

  "energy" : {
    "nonbonded" : {
      "coulomb" : { "epsr" : 1.0, "cutoff":9, "cutoffdip":4.3 },
      "ewald"   : { "precision":0.01, "epsr_surf":1e11, "cutoff":9, "alpha":0.2 },
      "cutoff_g2g" : 10.
    }
  },

  "atomlist" : {
    "OW" : { "q":-0.8476, "sigma":3.2, "eps":0.65 },
    "HW" : { "q": 0.4238, "sigma":0.0, "eps":0    },
    "dip" : { "sigma":1.0, "eps":0.17211, "mu":"0 0 0.415", "dp":0.3, "dprot":180,
              "alpha":"0 0 0 0 0 0", "theta":"0 0 0 0 0 0" }
  },

  "moleculelist" : {
    "sol" : { "atoms":"dip", "Ninit":100, "atomic":true },
    "water" : { "structure":"nemo.aam", "Ninit":216, "insdir":"1 1 1" }
    
  },

  "moves" : {
     "moltransrot" : {
        "water" : { "dp":0.5, "dprot":0.5, "dir":"1 1 1", "permol":true } 
     },
     "isobaric" : { "dp":0.1, "pressure":39.3155, "prop":1.0 },
    "atomtranslate" : { "sol" : { "peratom":true } },
    "atomrotate" : { "sol" : { "peratom":true } }
  },
  
  "system" : {
    "coulomb" : { "epsr":1, "cutoff":14., "cutoffdip":4.3 },
    "temperature"  : 300.0,
    "cuboid"   : { "len" : 27. },
    "mcloop"   : { "macro":10, "micro":100 },
    "atomlist"              : "nemo.json",
    "moleculelist"          : "nemo.json",
    "moleculecombinations"  : "nemo.json"
  }
  
}' > nemo.json

echo "3
OW    1    2.30    6.28    1.13 -0.8476  15.99940 1.6
HW    2    1.37    6.26    1.50  0.4238   1.00800 1
HW    3    2.31    5.89    0.21  0.4238   1.00800 1
" > nemo.aam


echo '{
  "atomlist" : {
    "dip" : { "sigma":1.0, "eps":0.17211, "mu":"0 0 0.415", "dp":0.3, "dprot":180,
              "alpha":"0 0 0 0 0 0", "theta":"0 0 0 0 0 0" }
  },

  "moleculelist" : {
    "sol" : { "atoms":"dip", "Ninit":100, "atomic":true }
  },

  "energy" : {
    "nonbonded" : {
       "coulomb" : { "epsr":1, "cutoff":4.3, "cutoffdip":4.3 }
    }
  },

  "moves" : {
    "atomtranslate" : { "sol" : { "peratom":true } },
    "atomrotate" : { "sol" : { "peratom":true } }
  },

  "system" : {
    "temperature"  : 300.,
    "cuboid"       : { "len" : 8.939 },
    "mcloop"       : { "macro":10, "micro":10000 },
    "atomlist"     : "nemo.json",
    "coulomb" : { "epsr":1, "cutoff":4.3, "cutoffdip":4.3 },
    "moleculelist" : "nemo.json"
  }
}' > nemo.json

echo '{
  "atomlist" : {
    "dip" : { "sigma":1.0, "eps":0.17211, "mu":"0 0 0.415", "dp":0.3, "dprot":180,
              "alpha":"0 0 0 0 0 0", "theta":"0 0 0 0 0 0" }
  },

  "moleculelist" : {
    "sol" : { "atoms":"dip", "Ninit":10000, "atomic":true }
  },

  "energy" : {
    "nonbonded" : {
       "coulomb" : { "epsr":1, "cutoff":4.3, "cutoff_rf":31.9165, "epsilon_rf":130.0 },
       "ewald"   : { "precision":0.01, "epsr_surf":1e11, "cutoff":19.0, "alpha":0.168421052631579, "maxK":11 },
    }
  },

  "moves" : {
    "atomtranslate" : { "sol" : { "peratom":true } },
    "atomrotate" : { "sol" : { "peratom":true } }
  },

  "system" : {
    "temperature"  : 300.,
    "cuboid"       : { "len" : 63.833 },
    "mcloop"       : { "macro":10, "micro":10 },
    "atomlist"     : "nemo.json",
    "coulomb" : { "epsr":1, "cutoff":4.3, "cutoffdip":4.3 },
    "moleculelist" : "nemo.json"
  }
}' > nemoRun.json


exe=./nemo
if [ -x $exe ]; then
 $exe
 rc=$?
 exit $rc
fi
exit 1
