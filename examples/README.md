# Faunus Examples

The `examples` directory contains a number of simple examples,
demonstrating features of Faunus. It is assumed that Faunus has
already been installed.

Example         | Description
--------------- | --------------------------------------------------------
`isochoric`     | Bending modulus via penalty function and isochoric moves
`membrane`      | Coarse-grained lipid bilayer under zero tension
`minimal`       | A simple Lennard-Jones particle system
`penalty`       | Wang-Landau sampling along a 2D reaction coordinate
`water`         | SPC/E water using Ewald or other electrostatic methods

Each example consists of a YAML input file and, potentially,
other files such as previously generated states (`.state.json`) or
output (`.out.json`). For instance, to run the `water` example, starting
from the stored state:

    $ yason.py water.yml | faunus --state water.state.json

This will run the simulation and generate the output file `out.json`
that can to converted to more readable, colorful YAML:

    $ yason.py out.json --color

