---
---
<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
});
</script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>
[![Edit](https://img.shields.io/badge/Github-Improve_this_page-orange.svg)]({{site.github.repository_url}}/blob/master/docs/{{page.path}})

# Monte Carlo Moves

A simulation can have an arbitrary number of MC moves operating on molecules, atoms, the volume, or
any other parameter affecting the system energy. Moves are specified in the `moves` section
at the top level input. For example:

~~~ yaml
moves:
    - moltransrot: { molecule: water, dp: 2.0, repeat: N
                    dprot: 1.0, dir: [1,1,0] }
    - volume: { dV: 0.01 }
    - ...

random: { seed: hardware }
~~~

The pseudo-random number engine used for MC moves can be seeded in three ways,

`seed`       | Description
-----------  | ----------------------------------------------
`fixed`      | Deterministic (default if `random` is absent)
`hardware`   | [Non-deterministric seed](http://en.cppreference.com/w/cpp/numeric/random/random_device)
engine state | [A previously saved stae](http://en.cppreference.com/w/cpp/numeric/random/mersenne_twister_engine/operator_ltltgtgt)

The last option is used to restore the state of the engine as saved along with normal simulation
output as a string containing a lenghty list of numbers.

## Translation and Rotation

### Molecular

`moltransrot`    |  Description
---------------- |  ---------------------------------
`molecule`       |  Molecule name to operate on
`dir=[1,1,1]`    |  Translational directions
`dp`             |  Translational displacement parameter
`dprot`          |  Rotational displacement parameter (radians)
`repeat=N`       |  Number of repeats per MC sweep. `N` equals $N_{molid}$ times.

This will simultaneously translate and rotate a molecular group by the following operation

$$
\textbf{r}^N_{trial} = \mbox{Rot}(\textbf{r}^N) + \delta
$$

where $\mbox{Rot}$ rotates `dprot`$\cdot \left (\zeta-\frac{1}{2} \right )$ radians around a random unit vector
emanating from the mass center,
$\zeta$ is a random number in the interval $[0,1[$, and
$\delta$ is a random unit vector scaled by `dp`.
Upon MC movement, the mean squared displacement
will be tracked.

### Atomic

`transrot`       |  Description
---------------- |  ---------------------------------
`molecule`       |  Molecule name to operate on
`dir=[1,1,1]`    |  Translational directions

As `moltransrot` but instead of operating on the molecular mass center, this translates
and rotates individual atoms in the group. The repeat is set to the number of atoms in the specified group and the displacement parameters `dp` and `dprot` for the individual atoms are taken from the atom properties defined in the [topology](../topology).

**note:**
atomic _rotation_ affects only anisotropic particles such as dipoles, spherocylinders, quadrupoles etc.
{: .notice--info}

### Cluster Move

`cluster`      | Description
-------------- | -----------------------
`molecules`    | List of molecules
`threshold`    | Mass-center threshold for forming a cluster
`dim=[1,1,1]`  | Directions to translate
`dprot`        | Rotational displacement
`dp`           | Translational displacement

This will attempt to rotate and translate clusters of molecular `molecules` defined by a distance `threshold`
between their mass centers. The move is associated with the following [bias](http://dx.doi.org/10/cj9gnn),
accounted for in the acceptance criterion:

## Pivot

`pivot`          | Description
---------------- | ----------------------------
`molecule`       | Molecule name to operate on
`dptot`          | Rotational displacement
`repeat=N`       | Number of repeats per MC sweep per bond

Performs a rotation around a random, harmonic bond vector in `molecule`, moving all atoms
either before _or_ after the bond with equal probability.

## Parallel Tempering

`temper`         | Description
---------------- | --------------------------------------------
`format=XYZQI`   | Particle properties to copy between replicas

We consider an extended ensemble, consisting of _n_
sub-systems or replicas, each in a distinct thermodynamic state (different
Hamiltonians) and with the total energy

$$
U = \sum_i^n\mathcal{H}_i(\mathcal{R}_i)
$$

The parallel tempering move performs a swap move where coordinate
spaces (positions, volume) between random, neighboring sub-systems, _i_ and _j_, are exchanged,

$$
\mathcal{R}_i^{\prime} = \mathcal{R}_j \quad \text{and} \quad \mathcal{R}_j^{\prime} = \mathcal{R}_i
$$

and the energy change of the _extended ensemble_, $\Delta U_{i\leftrightarrow j}$, is used in the
Metropolis acceptance criteria.

Parallel tempering requires compilation with MPI and the number
of replicas, _n_, exactly matches the number of processes. Each
replica prefixes input and output files with `mpi0.`, `mpi1.`,
etc. and only exchange between neighboring processes is performed.

**Note:**
Parallel tempering is currently limited to systems with
constant number of particles, $N$.
{: .notice--info}

## Volume Move <a name="volumemove"></a>

`volume`          |  Description
----------------- |  ----------------------------------------------
`dV`              |  Volume displacement parameter
`repeat=1`        |  Number of repeats per MC sweep.
`method=isotropic`|  Scaling method: `xy`, `isotropic`, `isochoric`

Performs a random walk in logarithmic volume,

$$
V^{\prime} = e^{\ln V + \left (\zeta-\frac{1}{2} \right )\cdot dV }
$$

and scales:

1. molecular mass centers
2. positions of free atoms (groups with `atomic=true`)

by $(V^{\prime}/V)^{1/3}$.
This is typically used for the $NPT$ ensemble, and for this an additional pressure term should be added to the Hamiltonian.
In the case of `isochoric` scaling, the total volume is kept constant and `dV` refers to an area change and reported output statistics on _volume_ should be regarded as _area_.
The table below explains the scaling behavior in different geometries:

`method`     |  Geometry    | Description
------------ | ------------ | ----------------------------
`isotropic`  |  `cuboid`    | Scales x, y, z
`isotropic`  |  `cylinder`  | Scales radius
`isotropic`  |  `sphere`    | Scales radius
`xy`         |  `cuboid`    | Scales xy, z untouched.
`isochoric`  |  `cuboid`    | Scales xy/z, const. volume


