[![Anaconda-Server Badge](https://anaconda.org/teokem/faunus/badges/installer/conda.svg)](https://conda.anaconda.org/teokem)
[![Anaconda-Server Badge](https://anaconda.org/teokem/faunus/badges/platforms.svg)](https://anaconda.org/teokem/faunus)[![Build Status](https://travis-ci.org/mlund/neofaunus.svg?branch=master)](https://travis-ci.org/mlund/neofaunus)

# Quick Start

- On mac or linux, install using Anaconda:

      $ conda install -c teokem faunus
    
  In addition to binaries in `bin`, examples are placed in `share/faunus`.
   
- To run the examples, the YAML input files must be piped through `yason.py`
  to create the required JSON input. For example,

      $ yason.py water.yml > water.json
      $ faunus -i water.json -o out.json

  Or simply:

      $ yason.py water.yml | faunus

  Reversely, JSON output can be converted to more readable YAML:

      $ yason.py out.json

- Documentation can be found [here](http://mlund.github.com/neofaunus/)
