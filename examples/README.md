# Faunus Examples

The `examples` directory contains a number of simple examples,
demonstrating features of Faunus. It is assumed that Faunus has
already been installed.

Example         | Description
--------------- | ---------------------------------------------------------
`bulk`          | Melted NaCl using Ewald summation and isobaric moves
`doublelayer`   | Electric double layer in a slit geometry with mean field correction
`cluster`       | Illustration of molecular cluster moves
`isochoric`     | Bending modulus via penalty function and isochoric moves
`membrane`      | Coarse-grained lipid bilayer under zero tension
`minimal`       | A simple Lennard-Jones particle system
`penalty`       | 2D Wang-Landau - see also Jupyter Notebook
`polymers`      | Hard-sphere polyelectrolytes in salt and NPT ensemble
`swapconf`      | Example of conformational swap moves from pool of weighted structures
`water`         | SPC/E water in the NPT ensemble

Each example consists of a YAML input file and, potentially,
other files such as previously generated states (`.state.json`) or
output (`.out.json`). For instance, to run the `water` example, starting
from the stored state:

    $ yason.py water.yml | faunus --state water.state.json

This will run the simulation and generate the output file `out.json`
that can to converted to more readable, colorful YAML:

    $ yason.py out.json --color

