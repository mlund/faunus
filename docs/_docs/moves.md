<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
});
</script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

# Monte Carlo Moves

A simulation can have an arbitrary number of MC moves operating on molecules, atoms, the volume, or
any other parameter affecting the system energy. Moves are specified in the `moves` section
at the top level input. For example:

~~~ yaml
moves:
    - moltransrot: { molecule: water, dp: 2.0, repeat: N,
                    dprot: 1.0, dir: [1,1,0] }
    - volume: { dV: 0.01 }
    - ...

random: { seed: hardware }
~~~

The pseudo-random number engine used for MC moves can be seeded in three ways,

`seed`       | Description
------------ | ----------------------------------------------
`fixed`      | Deterministic (default if `random` is absent)
`hardware`   | [Non-deterministric seed](http://en.cppreference.com/w/cpp/numeric/random/random_device)
engine state | [A previously saved stae](http://en.cppreference.com/w/cpp/numeric/random/mersenne_twister_engine/operator_ltltgtgt)

The last option is used to restore the state of the engine as saved along with normal simulation
output as a string containing a lenghty list of numbers.
If initialization from a previously saved state fails -- this may happen if
generated on another operating system -- a warning is issued and the seed
falls back to `fixed`.

## Translation and Rotation

The following moves are for translation and rotation of atoms, molecules, or clusters.
The `dir` keyword restricts translational directions which by default is set to
`[1,1,1]`, meaning translation by a unit vector, randomly picked on a sphere, and scaled by
a random number in the interval `[0, dp]`.
If `dir=[1,1,0]` the unit vector is instead picked on a circle (here _x_, _y_)
and if `dir=[0,0,1]` on a line (here $z$).

### Molecular

`moltransrot`    |  Description
---------------- |  ---------------------------------
`molecule`       |  Molecule name to operate on
`dir=[1,1,1]`    |  Translational directions
`dirrot=[0,0,0]` |  Predefined axis of rotation
`dp`             |  Translational displacement parameter
`dprot`          |  Rotational displacement parameter (radians)
`repeat=N`       |  Number of repeats per MC sweep. `N` equals $N\_{molid}$ times.

This will simultaneously translate and rotate a molecular group by the following operation

$$
\textbf{r}^N\_{trial} = \mbox{Rot}(\textbf{r}^N) + \delta
$$

where $\mbox{Rot}$ rotates `dprot`$\cdot \left (\zeta-\frac{1}{2} \right )$ radians around a random unit vector
emanating from the mass center,
$\zeta$ is a random number in the interval $[0,1[$, and
$\delta$ is a random unit vector scaled by a random number in the interval `[0,dp]`.
A predefined axis of rotation can be specified as `dirrot`. For example, setting `dirrot` to [1,0,0], [0,1,0] or [0,0,1] 
results in rotations about the $x-$, $y-$, and $z-$axis, respectively.
Upon MC movement, the mean squared displacement will be tracked.


### Atomic

`transrot`          |  Description
------------------- |  ---------------------------------
`molecule`          |  Molecule name to operate on
`dir=[1,1,1]`       |  Translational directions
`energy_resolution` |  If set to a non-zero value (kT), an energy histogram will be generated.

As `moltransrot` but instead of operating on the molecular mass center, this translates
and rotates individual atoms in the group. The repeat is set to the number of atoms in the specified group and the
displacement parameters `dp` and `dprot` for the individual atoms are taken from
the atom properties defined in the [topology](topology).
Atomic _rotation_ affects only anisotropic particles such as dipoles, spherocylinders, quadrupoles etc.

An energy histogram of each participating species will be written to disk if the `energy_resolution`
keyword is set. The value (in kT) specifies the resolution of the histogram binning. The analysis is
essentially for free as the energies are already known from the move.


### Cluster Move

`cluster`       | Description
--------------- | --------------------------------------------
`molecules`     | Array of molecule names; `[*]` selects all
`threshold`     | Mass-center threshold for forming a cluster (number or object)
`dir=[1,1,1]`   | Directions to translate
`dirrot=[0,0,0]`| Predefined axis of rotation. If zero, a random unit vector is generated for each move event
`dprot`         | Rotational displacement (radians)
`dp`            | Translational displacement (Å)
`single_layer=false` | If `true`, stop cluster-growth after one layer around centered molecule (experimental)
`satellites`    | Subset of `molecules` that cannot be cluster centers
`com=true`      | Use distance threshold between mass-centers instead of particles when finding clusters
`com_shape=true`| Use mass centers for shape analysis instead of particle positions (affects analysis only)
`analysis`      | See below

This will attempt to rotate and translate clusters of molecular `molecules` defined by a distance `threshold`
between mass centers.
The `threshold` can be specified as a single distance or as a complete list of combinations, see example below.
For simulations where small molecules cluster around large macro-molecules, it can be useful to use the `satellites`
keyword which denotes a list of molecules that can be part of a cluster, but cannot be the cluster nucleus or
starting point.
All molecules listed in `satellites` must be part of `molecules`.
A predefined axis of rotation can be specified as `dirrot`. For example, setting `dirrot` to [1,0,0], [0,1,0] or [0,0,1] 
results in rotations about the $x-$, $y-$, or $z-$axis, respectively.

The move is associated with [bias](http://dx.doi.org/10/cj9gnn), such that
the cluster size and composition remain unaltered.
If a cluster is larger than half the simulation box length, only translation will be attempted.

Example:

``` yaml
cluster:
   molecules: [protein, cations]
   satellites: [cations]
   threshold:
      protein protein: 25
      protein cations: 15
      cations cations: 0
   dp: 3
   dprot: 1
   analysis: {file: cluster_shape.dat.gz}
```

Using `analysis` the move also reports on the average cluster size; cluster size distribution; and
relative shape anisotropy, $\kappa^2$ as a function of cluster size. If all particles in the cluster
are distributed on a sphere then $\kappa^2=0$, while if on a straight line, $\kappa^2=1$. See the
`polymershape` analysis for more information.

`analysis`    | Description
------------- | ----------------------------------------------------------------------------
`com=true`    | Use molecule mass center instead of particle positions for shape anisotropy
`file`        | If given save shape properties for each sample event
`save_pqr`    | If set to true, save PQR files of cluster, one for each size 
`interval=10` | Interval between samples


### Smarter Monte Carlo

Preferential selction of particles can be enabled via the `region` keyword which instructs
some moves to pick particles or groups preferentially from a given _region_. As described
in [doi:10/frvx8j](https://doi.org/frvx8j) a bias is introduced which is automatically
accounted for. The preference for sampling inside the region is controlled by `p` which
can be regarded as an outside update probability.
If $p=1$ no preferential sampling is performed, whereas if
$p<1$, more sampling is done inside the region.

For example:

~~~ yaml
- moltransrot:
    ...
    ...
    region:
      policy: ellipsoid
      parallel_radius: 5.0
      perpendicular_radius: 4.0
      index1: 10
      index2: 12
      p: 0.2
~~~

The available regions are:

#### Ellipsoid

The connection vector between two (moving) reference particles defines an ellipsoid
centered at the midpoint between the reference particles.
The reference particle separation is unimportant, only the direction is used.

`policy=ellipsoid`     | Description
---------------------- | ----------------------------------------------------------------
`p`                    | Number (0,1] where a lower number means higher regional sampling
`index1`               | Index of first reference particle
`index2`               | Index of second reference particle
`parallel_radius`      | Radius parallel to axis connecting the two references
`perpendicular_radius` | Radius perpendicular to axis connecting the two references
`group_com=false`      | Test only mass center of molecular groups


#### Within Molecule (experimental)

Samples from within a threshold from a molecule type. This can be useful to _e.g._
preferentially sample solvent molecules around dilute solute molecules.
The `com` keyword is available if the selected `molecule` has a well-defined mass-center,
_i.e._ if `is_atomic=false`.
It is also possible to use only the mass center for the moved groups by setting `group_com`.

`policy=within_molid`  | Description
---------------------- | ----------------------------------------------------------------
`p`                    | Number (0,1] where a lower number means higher regional sampling
`molecule`             | Name of molecule to search around
`threshold`            | Distance threshold to any particle in `molecule`
`com=false`            | Use `threshold` with respect to mass-center of `molecule`
`group_com=false`      | Test with respect to mass center of molecular groups


## Internal Degrees of Freedom

### Charge Move

`charge`         | Description
---------------- | ---------------------------------
`index`          | Atom index to operate on
`dq`             | Charge displacement
`quadratic=true` | Displace linearly along q^2 instead of q

This performs a fractional charge move on a specific atom. The charge
displacement can be performed linerly along $q$ or linearly along $q^2$.
For the latter the following bias energy will be added to ensure
uniform sampling of $q$,

$$
u = k\_BT\ln \left ( \left | q^{\prime} / q \right |\right )
$$

Limitations:
This move changes the particle charge and therefore cannot be used with
splined pair-potentials where the initial charges from are read from `atomlist`.
Instead, use a hard-coded variant like `nonbonded_coulomblj` etc.


### Conformational Swap

`conformationswap`  | Description
------------------- | ---------------------------------
`molecule`          | Molecule name to operate on
`repeat=N`          | Number of repeats per MC sweep
`keeppos=False`     | Keep original positions of `traj`
`copy_policy=all`   | What to copy from library. See table below.

This will swap between different molecular conformations
as defined in the [Molecule Properties](topology.html#molecule-properties) with `traj` and `trajweight`
If defined, the weight
distribution is respected, otherwise all conformations
have equal intrinsic weight. Upon insertion, the new conformation
is randomly oriented and placed on top of the mass-center of
an exising molecule. That is, there is no mass center movement.
If `keeppos` is activated, the raw coordinates from the conformation
is used, i.e. no rotation and no mass-center overlay.

By default all information from the conformation is copied (`copy_policy=all`), including charges and particle type.
To for example copy only positions, use `copy_policy=positions`. This can be useful when using speciation moves.

`copy_policy`  | What is copied
-------------- | ----------------------------------------------------
`all`          | All particle properties
`positions`    | Positions, only
`charges`      | Charges, only
`patches`      | Spherocylinder patch and length, but keep directions


### Pivot

`pivot`          | Description
---------------- | ----------------------------
`molecule`       | Molecule name to operate on
`dprot`          | Rotational displacement
`repeat=N`       | Number of repeats per MC sweep
`skiplarge=true` | Skip too large molecules

Performs a rotation around a random, harmonic bond vector in `molecule`, moving all atoms
either before _or_ after the bond with equal probability. Current implementation assumes
unbranched chains with all atoms as links, i.e., no side chains are present.
For long polymers (compared to the box size), a large displacement parameter may cause
problems with mass center calculation in periodic systems.
This can be caught with the `sanity` analysis and should it occur, try one of the following:

- enable `skiplarge`
- decrease `dprot`
- increase the simulation container.

The first option will simply reject troublesome configurations and the
final output contains information of the skipped fraction. Skipping is
unphysical so make sure the skipped fraction is small.

The default value of `repeat` is the number of harmonic bonds in the `molecule`
(multiplied by the number of molecules).

Limitations: Chain bonds have to be ordered sequentially in the topology.


### Crankshaft

`crankshaft`         | Description
-------------------- | --------------------------------------------------------
`molecule`           | Molecule name to operate on
`dprot`              | Rotational displacement
`repeat=N`           | Number of repeats per MC sweep
`skiplarge=true`     | Skip too large molecules
`joint_max=`$\infty$ | Maximum number of bonds between randomly selected atoms


Performs a rotation of a chain segment between two randomly selected atoms in the `molecule`.

The default value of `repeat` is the number of atoms in the `molecule` minus two
(multiplied by the number of molecules).


## Parallel Tempering

`temper`                 | Description
------------------------ | ----------------------------------------------------------------------
`format=xyzqi`           | Particle properties to copy between replicas (`xyzqi`, `xyzq`, `xyz`)
`volume_scale=isotropic` | How to apply exchanged volumes: `z`, `xy`, `isotropic`, `isochoric`
`nstep=1`                | Number of sweeps between samples.
`partner_policy=oddeven` | Policy used to create partner pairs (currently only `oddeven`)

We consider an extended ensemble, consisting of _n_
sub-systems or replicas, each in a distinct thermodynamic state (different
Hamiltonians) and with the total energy

$$
U = \sum\_i^n\mathcal{H}\_i(\mathcal{R}\_i)
$$

The parallel tempering move performs a swap move where coordinate
spaces (positions, volume) between random, neighboring sub-systems, _i_ and _j_, are exchanged,

$$
\mathcal{R}\_i^{\prime} = \mathcal{R}\_j \quad \text{and} \quad \mathcal{R}\_j^{\prime} = \mathcal{R}\_i
$$

and the energy change of the _extended ensemble_, $\Delta U\_{i\leftrightarrow j}$, is used in the
Metropolis acceptance criteria.

Parallel tempering requires compilation with MPI and the number
of replicas, _n_, exactly matches the number of processes. Each
replica prefixes input and output files with `mpi0.`, `mpi1.`,
etc. and only exchange between neighboring processes is performed.
The move is is performed exactly every `nstep` Monte Carlo cycle.
By default, particle positions (`xyz`), charge (`q`), and atom id (`i`) are exchanged
between replicas and can be controlled with `format`.

Support for fluctuating number of particles, i.e.
grand canonical moves is currently untested and should be
considered experimental.

#Exchange statistics
In the output file, the acceptance and attempts of an exchange direction can be found under exchange in temper;moves. Next, the first tap under  exchangempi gives which mpi box the outfile represents (ex. 0, 1,...). Under this number the tab exchanges all executed temper movements are listed. For each movement two placeholders  are given, which represents the mpi's the tempering algorithm tries to switch. If the tempering move was not accepted, both the placeholders will be -1. When the tempering move is accepted do the placeholders contain the numbers of the interchanged mpi's. In faunus/examples/temper there is a notebook "temper\_exchange\_statistics.ipynb" illustrating how this aquired data can be used to study how the mpi's are exchanged.  


## Volume Move

`volume`          |  Description
----------------- |  ----------------------------------------------
`dV`              |  Volume displacement parameter
`repeat=1`        |  Number of repeats per MC sweep.
`method=isotropic`|  Scaling method: `z`, `xy`, `isotropic`, `isochoric`

Performs a random walk in logarithmic volume,

$$
V^{\prime} = e^{\ln V + \left (\zeta-\frac{1}{2} \right )\cdot dV }
$$

and scales:

1. molecular mass centers
2. positions of free atoms (groups with `atomic=true`)

by $(V^{\prime}/V)^{1/3}$.
This is typically used for the $NPT$ ensemble, and for this an additional pressure term should be added to the Hamiltonian.
In the case of `isochoric` scaling, the total volume is kept constant and `dV` refers to an area change and reported output
statistics on _volume_ should be regarded as _area_.
The table below explains the scaling behavior in different geometries:

`method`     |  Geometry    | Description
------------ | ------------ | ----------------------------
`isotropic`  |  `cuboid`    | Scales x, y, z
`isotropic`  |  `cylinder`  | Scales radius
`isotropic`  |  `sphere`    | Scales radius
`z`          |  `cuboid`    | Scales z, xy untouched.
`xy`         |  `cuboid`    | Scales xy, z untouched.
`isochoric`  |  `cuboid`    | Scales xy/z, const. volume

For cuboidal geometries, the scaling in each of the specified dimensions is $(V^{\prime}/V)^{1/d}$,
where $d=3$ for `isotropic`, $d=2$ for `xy`, and $d=1$ for `z`.

_Warning:_ Untested for cylinders, slits.

## Gibbs Ensemble (unstable)

_Note: this is marked unstable or experimental, meaning that it is still being tested and
may see significant changes over time._

[Gibbs ensemble](https://dx.doi.org/10/cvzgw9)
can be used to investigate phase transitions by _matter_ and _volume_ exchange between two cells.
The `examples/gibbs-ensemble/` directory contains a Jupyter Notebook with a worked example of a simple Lennard-Jones system.
Multi-component mixtures are supported via the required `molecules` and `molecule` keywords which indicate which species
are to be affected.
Volume and matter exchange are done in separate moves, the latter _per_ species:

~~~ yaml
insertmolecule:
  - A: 100, inactive: 50} # note inactive species
  - B: 100, inactive: 50}
moves:
  - gibbs_volume: { dV: 1.0, molecules: ["A", "B"] } # exchange volume
  - gibbs_matter: { molecule: "A" } # exchange A molecules
  - gibbs_matter: { molecule: "B" } # exchange B molecules
~~~

In addition, you will likely also want to add translational and rotational moves.
It is important that each cell can accommodate _all_ particles in the system.
This is done by reserving an appropriate number of `inactive` particles in the initial
configuration, see above example.
An error is thrown if this criterion is neglected.

### Running

Gibbs ensemble requires that Faunus is compiled with MPI support, check with `faunus --version`,
and _exactly two_ processes must be give with e.g. `mpirun -np 2`.

- If starting conditions for each cell are identical, use `--nopfx` and a single `input.json` file:
  ~~~ bash
  mpirun -np 2 faunus --nopfx --input input.json
  ~~~
- If input differs, e.g. different initial volumes or number of particles, create two input files, prefixed with `mpi0.` and `mpi1.`,
  and skip the `--nopfx` flag.
- Reload from existing states by using the `--state` flag. `mpi` prefix are automatically added.


## Reactive Canonical Monte Carlo

The speciation move handles density fluctuations and particle transformations and is the main move
for particle insertion, deletion, and swapping used in (semi)-grand canonical ensembles.
A reaction from `reactionlist` is randomly picked from the topology and is either propagated forward or backward.
In Faunus, the total number of atoms and molecules is constant, but these can be either _active_ or _inactive_.
Deleting a molecule simply deactivates it, while insertion _vice versa_ activates an inactive molecule.
Thus, it is important that the _capacity_ or reservoir of particles (active plus inactive) is
sufficiently large to allow for fluctuations.
This is ensured using `insertmolecules` (see Topology).
A runtime warning will be given, should you run low on particles.<br>
Besides deleting/inserting molecules (mono- or polyatomic), the speciation move performs reactions involving a
single-atom ID transformation (_e.g._, acid-base reactions).
In this case, an particle of type A (part of a mono- or polyatomic molecule) is randomly picked from the system
and all its properties, except its position, are replaced with those of an atom of type B.
Such ID transormations can also involve the addition/deletion of molecules or _implicit_ atoms.<br>
For a reaction
$$
\sum\_i \nu\_i M\_i = 0
$$
where $M\_i$ is the chemical symbol and $\nu\_i$ is the stoichiometric coefficient of species $i$ (positive for products and negative for reagents),
the contribution of a speciation move to the energy change is
$$
\beta \Delta U = - \sum\_i \ln{ \left ( \frac{ N\_i! }{(N\_i+\nu\_i)!} V^{\nu\_i} \right ) } - \ln{ \prod\_i a\_i^{\nu\_i} },
$$
where $N\_i$ is the number of particles of species $i$ in the current state and $a\_i$ is the activity of species $i$.

For more information, see the Topology section and [doi:10/fqcpg3](https://doi.org/10/fqcpg3).

`rcmc`          |  Description
--------------- | ----------------------------------
`repeat=1`      |  Average number of moves per sweep


## Replay

`replay`         | Description
---------------- | ----------------------------
`file`           | Trajectory file to read (xtc)

Use next frame of the recorded trajectory as a move. The move is always unconditionally accepted,
hence it may be used to replay a simulation, e.g., for analysis. Currently only Gromacs compressed
trajectory file format (XTC) is supported. Note that total number of steps (macro × micro) should
correspond to the number of frames in the trajectory.
