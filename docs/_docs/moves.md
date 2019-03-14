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
$\delta$ is a random unit vector scaled by a random number in the interval `[0,dp]`.
Upon MC movement, the mean squared displacement
will be tracked.


### Atomic

`transrot`       |  Description
---------------- |  ---------------------------------
`molecule`       |  Molecule name to operate on
`dir=[1,1,1]`    |  Translational directions

As `moltransrot` but instead of operating on the molecular mass center, this translates
and rotates individual atoms in the group. The repeat is set to the number of atoms in the specified group and the
displacement parameters `dp` and `dprot` for the individual atoms are taken from
the atom properties defined in the [topology](../topology).

**note:**
atomic _rotation_ affects only anisotropic particles such as dipoles, spherocylinders, quadrupoles etc.
{: .notice--info}

### Cluster Move

`cluster`      | Description
-------------- | -----------------------
`molecules`    | Array of molecule names; `[*]` selects all 
`threshold`    | Mass-center threshold for forming a cluster
`dir=[1,1,1]`  | Directions to translate
`dprot`        | Rotational displacement (radians)
`dp`           | Translational displacement

This will attempt to rotate and translate clusters of molecular `molecules` defined by a distance `threshold`
between their mass centers. The move is associated with [bias](http://dx.doi.org/10/cj9gnn), such that
the cluster size and composition remain unaltered.
If a cluster is larger than half the simulation box length, only translation will be attempted.

**Restrictions:**
Currently, the number of `molecules` must be constant throughout simulation.
{: .notice--info}



## Internal Degrees of Freedom

### Charge Move

`charge`         |  Description
---------------- |  ---------------------------------
`index`          |  Atom index to operate on
`dq`             |  Charge displacement

This performs a fractional charge move on a specific atom.


### Conformational Swap

`conformationswap` | Description
------------------ | ---------------------------------
`molecule`         |  Molecule name to operate on
`repeat=N`         |  Number of repeats per MC sweep

This will swap between different molecular conformations
as defined in the topology with `traj` and `trajweight`
If defined, the weight
distribution is respected, otherwise all conformations
have equal intrinsic weight. Upon insertion, the new conformation
is randomly oriented and placed on top of the mass-center of
an exising molecule. That is, there is no mass center movement.


### Pivot

`pivot`          | Description
---------------- | ----------------------------
`molecule`       | Molecule name to operate on
`dprot`          | Rotational displacement
`repeat=N`       | Number of repeats per MC sweep per bond
`skiplarge=true` | Skip too large molecules

Performs a rotation around a random, harmonic bond vector in `molecule`, moving all atoms
either before _or_ after the bond with equal probability.
For long polymers (compared to the box size), a large displacement parameter may cause
problems with mass center calculation in periodic systems.
This can be caught with the `sanity` analysis and should it occur, try one of the following:

- enable `skiplarge`
- decrease `dprot`
- increase the simulation container.

The first option will simply reject troublesome configurations and the
final output contains information of the skipped fraction. Skipping is
unphysical so make sure the skipped fraction is small.


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

**Warning:**
Untested for cylinders, slits.
{: .notice--warning}


## Reactive Canonical Monte Carlo (beta) <a name="rcmc"></a>

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
\sum_i \nu_i M_i = 0
$$
where $M_i$ is the chemical symbol and $\nu_i$ is the stoichiometric coefficient of species $i$ (positive for products and negative for reagents),
the contribution of a speciation move to the energy change is
$$
\beta \Delta U = - \sum_i \ln{ \left ( \frac{ N_i! }{(N_i+\nu_i)!} V^{\nu_i} \right ) } - \ln{ \prod_i a_i^{\nu_i} },
$$
where $N_i$ is the number of particles of species $i$ in the current state and $a_i$ is the activity of species $i$.

For more information, see the Topology section and [doi:10/fqcpg3](https://doi.org/10/fqcpg3).

`rcmc`          |  Description
--------------- | ----------------------------------
`repeat=1`      |  Average number of moves per sweep

**Warning:**
The speciation move is under construction and subject to change.
{: .notice--warning}

