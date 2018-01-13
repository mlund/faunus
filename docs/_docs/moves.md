---
---
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
    - moltransrot: { molecule: water, dp: 2.0, repeat: N
                    dprot: 1.0, dir: [1,1,0] }
    - volume: { dV: 0.01 }
    - ...

random:
    seed: hardware
~~~

The pseudo-random number engine used for MC moves can be seeded in three ways,

`seed`       | Description
-----------  | ----------------------------------------------
`default`    | Deterministic (default if `random` is absent)
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

**Note:**
Atomic _rotation_ affects only anisotropic particles such as dipoles, spherocylinders, quadrupoles etc.
{: .notice--info}

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

Here we consider an extended ensemble, consisting of a _n_
sub-systems (replicas), each in a different thermodynamic state (different
Hamiltonians) and with a total energy

$$
U = \sum_i^n\mathcal{H}_i(\mathcal{R_i})
$$

The parallel tempering move performs a swap move where the coordinate
space, $\mathcal{R}$, between two sub-systems is exchanged and
$\Delta U_{i\leftrightarrow j}$ is the resulting trial energy.

**Note:**
Parallel tempering require compilation with MPI and the number
of replicas, _n_, exactly match the number of processes. Each
replica will prefix input and output files with `mpi0.`, `mpi1.`,
etc.
{: .notice--info}

## Volume Move <a name="volumemove"></a>

`volume`         |  Description
---------------- |  ---------------------------------
`dV`             |  Volume displacement parameter
`repeat=1`       |  Number of repeats per MC sweep.

Performs a random walk in logarithmic volume,

$$
V^{\prime} = e^{\ln V + \left (\zeta-\frac{1}{2} \right )\cdot dV }
$$

and scales:

1. molecular mass centers
2. positions of free atoms (groups with `atomic=true`)

by $(V^{\prime}/V)^{1/3}$. This is typically used for the $NPT$ ensemble, and for this an additional pressure term should be added to the Hamiltonian.
