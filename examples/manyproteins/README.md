# Description of system

Simulation of 300 rigid, molecules (globular proteins) in a cubic box.

- each molecule has 272 interaction sites
- sampling is done using Monte Carlo, where ONE protein is moved at the time, followed by an energy difference calculation (bottleneck)
- the molecule-molecule interaction is a N^2 operation
- the site interaction potential is constructed from two functions (`energy -> nonbonded`) that are added together and SPLINED.
  see `Potential::SplinedPotential`.
- a mass-center cutoff scheme is applied which mimicks a cell-list
- use `energy -> summation_policy` to control (parallel) energy summation schemes

## Code information

The software, `faunus`, is programmed in C++17 and all project files are in `src/`. Supporting libraries
are downloading and build during the `cmake` configure step. 
The Matrix library `Eigen` is used for all vector operations, but we mainly use cartesian distances.

## Suspected bottle necks / targets

- Summing of pair-wise interactions between sites (N^2 complexity).
- `tabulate.h` --> `Andrea::eval()` function for splining (SIMD optimization?)
- `energy.h` --> `DelayedEnergyAccumulator` performs the actual summing; already contains OMP pragmas / or STL parallel algorithms.
- `geometry.h` --> `Cuboid::sqdist` and `Cuboid::vdist` for distance calculations (SIMD optimization?)
- `move.h` --> `perform_move()`. Here the old and new energy may be computed in two `omp sections`, but seem for hamper contained
  `omp parallel for`.

## Running the example

Requires that `ruamel_yaml` and `jsonschema` are installed, typically using conda.
The `yason.py` script simply converts from YAML to JSON and at the same time checks the input syntax.

~~~ bash
cd examples/manyproteins
../../faunus --input input.json --verbosity 5 --state restart.json 
~~~
