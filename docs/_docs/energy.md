---
---
<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
});
</script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>
[![](https://img.shields.io/badge/Github-Improve_this_page-orange.svg)]({{site.github.repository_url}}/blob/master/docs/{{page.path}})

# Energy <a name="energy"></a>

The system energy, or Hamiltonian, consists of a sum of potential energy terms,

$$
\mathcal{H}_{sys} = U_1 + U_2 + ...
$$

The energy terms are specified in `energy` at the top level input.
For example:

~~~ yaml
energy:
    - isobaric: {P/atm: 1}
    - sasa: {molarity: 0.2, radius: 1.4 }
    - confine: {type: sphere, radius: 10,
                 molecules: [water]}
    - nonbonded:
        default: # applied to all atoms
        - lennardjones: {mixing: LB}
        - coulomb: {type: plain, epsr: 1, cutoff: 12}
        Na CH:   # overwrite specific atom pairs
        - wca: { mixing: LB }

    - ...
~~~

**Note:**
_Energies_ in MC may contain implicit degrees of freedom, _i.e._ be temperature-dependent,
effective potentials. This is inconsequential for sampling
density of states, but care should be taken when interpreting derived functions such as
energies, entropies, pressure etc.
{: .notice--info}

## External Pressure <a name="isobaric"></a>

This adds the following pressure term[^frenkel] to the Hamiltonian, appropriate for
MC moves in $\ln V$:

$$
    U = PV - k_BT\left ( N + 1 \right ) \ln V
$$

where $N$ is the total number of molecules and atomic species.

[^frenkel]: _Frenkel and Smith, 2nd Ed., Chapter 5.4_.

`isobaric`   | Description
-------------|-------------------------------------------------------
`P/unit`     | External pressure where unit can be `mM`, `atm`, or `Pa`.


## Nonbonded Interactions <a name="nonbonded"></a>

This term loops over pairs of atoms, $i$, and $j$, summing a given pair-wise additive potential, $u_{ij}$,

$$ U = \sum_{i=0}^{N-1}\sum_{j=i+1}^N u_{ij}(\textbf{r}_j-\textbf{r}_i)$$

**Note:** the pair-potential can be a combination of several potentials defined at runtime.
However, for optimal performance we include certain hard-coded combinations, defined at _compile time_.
It is straight forward to add more by editing the class
`Energy::Hamiltonian` found in `src/energy.h` and then re-compile.
{: .notice--info}

Below is a description of possible pair-potentials and their configuration.

`energy`               | $u_{ij}$
---------------------- | ----------------------------------
`nonbonded`            | Any combination of pair potentials
`nonbonded_coulomblj`  | `coulomb`+`lennardjones`
`nonbonded_coulombhs`  | `coulomb`+`hardsphere`
`nonbonded_coulombwca` | `coulomb`+`wca`
`nonbonded_pmwca`      | `coulomb`+`wca` (`type=plain`, `cutoff`$=\infty$)


### Electrostatics

 `coulomb`   |  Description
 ----------- |  -------------------------------------------------
 `type`      |  Coulomb type, see below
 `cutoff`    |  Spherical cutoff, $R_c$ after which the potential is zero
 `epsr`      |  Relative dielectric constant of the medium
 `utol=1e-5` |  Error tolerence for splining

This is a multipurpose potential that handles several electrostatic methods.
Beyond a spherical real-space cutoff, $R_c$, the potential is zero while if
below,

$$
u_{ij} = \frac{e^2 z_i z_j }{ 4\pi\epsilon_0\epsilon_r r_{ij} }\mathcal{S}(q)
$$

where $\mathcal{S}(q=r/R_c)$ is a splitting function:

`type`                                   | Keywords      | $\mathcal{S}(q)$                                   
---------------------------------------- | ------------- | ---------------------------------------------------
`none`                                   |               | 0
[`plain`](http://doi.org/ctnnsj)         |               | 1
[`fanourgakis`](http://doi.org/f639q5)   |               | $1-\frac{7}{4}q+\frac{21}{4}q^5-7q^6+\frac{5}{2}q^7$
[`ewald`](http://doi.org/dgpdmc)         | `alpha`       | $\text{erfc}(\alpha R_cq)$
[`wolf`](http://doi.org/cfcxdk)          | `alpha`       | $\text{erfc}(\alpha R_cq)-\text{erfc}(\alpha R_c)q$
[`yukawa`](http://bit.ly/2CbVJ3v)        | `debyelength` | $e^{-\kappa R_c q}-e^{-\kappa R_c}$
[`yonezawa`](http://dx.doi.org/10/j97)   | `alpha`       | $1+\text{erfc}(\alpha R_c)q+q^2$
[`qpotential`](http://goo.gl/hynRTS)     | `order=300`   | $\prod_{n=1}^{\text{order}}(1-q^n)$
[`reactionfield`](http://doi.org/dbs99w) | `epsrf`       | $1+\frac{\epsilon_{RF}-\epsilon_r}{2\epsilon_{RF}+\epsilon_r}q^3-3\frac{\epsilon_{RF}}{2\epsilon_{RF}+\epsilon_r}q$
[`fennel`](http://doi.org/bqgmv2)        | `alpha`       | $\scriptstyle\text{erfc}(\alpha R_cq)-\text{erfc}(\alpha R_c)q+(q-1)q \left( \text{erfc}(\alpha R_c) + \frac{2\alpha R_c}{\sqrt{\pi}} e^{-\alpha^2 R_c^2} \right)$

**Note:** Internally $\mathcal{S}(q)$ is _splined_ whereby all types evaluate at similar speed.
{: .notice--info}

#### Ewald Summation

If type is `ewald`, terms from reciprocal space; surface energies; and
self energies are automatically added to the Hamiltonian, activating additional keywords:

`type=ewald`         | Description
-------------------- | ---------------------------------------------------------------------
`kcutoff`            | Reciprocal-space cutoff
`epss=0`             | Dielectric constant of surroundings, $\varepsilon_{surf}$ (0=tinfoil)
`ipbc=false`         | Use isotropic periodic boundary conditions, IPBC.
`spherical_sum=true` | Spherical/ellipsoidal summation in reciprocal space; cubic if `false`.

The added energy terms are:

$$
\small
\begin{aligned}
U =& \overbrace{\frac{2\pi f}{V}\sum_{ {\bf k} \ne {\bf 0}} A_k\vert Q^{q\mu} \vert^2}^{\text{reciprocal}}
- \overbrace{ f \sum_{j} \left( \frac{\alpha}{\sqrt{\pi}}q_j^2 + \frac{2\alpha^3}{3\sqrt{\pi}}\vert{\boldsymbol{\mu}}_j\vert^2   \right)}^{\text{self}}\\
&+ \underbrace{\frac{2\pi f}{(2\varepsilon_{surf} + 1)V}\left(  \vert \sum_{j}q_j{\bf r}_j   \vert^2 + 2\sum_{j}q_i{\bf r}_j \cdot \sum_{j}{\boldsymbol{\mu}}_j + \vert \sum_{j}{\boldsymbol{\mu}}_j \vert^2 \right )}_{\text{surface}}\\
\end{aligned}
$$

where

$$
    f = \frac{1}{4\pi\varepsilon_0\varepsilon_r} \quad\quad V=L_xL_yL_z
$$

$$
A_k = \frac{e^{-k^2/4\alpha^2}}{k^2}
\quad \quad Q^{q\mu} = \sum_{j}q_j + i({\boldsymbol{\mu}}_j\cdot {\bf k})  e^{i({\bf k}\cdot {\bf r}_j)}
$$

$$
{\bf k} = 2\pi\left( \frac{n_x}{L_x} , \frac{n_y}{L_y} ,\frac{n_z}{L_z} \right),\;\; {\bf n} \in \mathbb{Z}^3
$$

In the case of isotropic periodic boundaries (`ipbc=true`), the orientational degeneracy of the
periodic unit cell is exploited to mimic an isotropic environment, reducing the number
of wave-vectors by one fourth compared with PBC Ewald.
For point charges, IPBC introduce the modification,

$$
Q^q = \sum_jq_j\prod_{\alpha\in\{x,y,z\}}\cos\left(\frac{2\pi}{L_{\alpha}}n_{\alpha} r_{\alpha,j}\right)
$$

while for point dipoles,

$$
Q^{\mu} = \sum_j\boldsymbol{\mu}_j\cdot\nabla_j\left(\prod_{\alpha \in\{x,y,z\}}\cos\left(\frac{2\pi}{L_{\alpha}}n_{\alpha}r_{\alpha,j}\right)\right).
$$

**Limitations:** Ewald summation requires a constant number of particles, i.e. $\mu V T$ ensembles
and Widom insertion are currently unsupported.
{: .notice--info}

### Charge-Nonpolar

The energy when the field from a point charge, $z_i$, induces a dipole in a polarizable particle of unit-less excess polarizability, $\alpha_j=\left ( \frac{\epsilon_j-\epsilon_r}{\epsilon_r+2\epsilon_r}\right ) a^3$, is

$$
    \beta u_{ij} = -\frac{\lambda_B z_i^2 \alpha_j}{2r_{ij}^4}
$$

where $a_j$ is the radius of the non-polar particle and $\alpha_j$ is set in
the atom topology, `alphax`.
If `alphaneutral=true` it is required that one of the particles
is charged, while the other is neutral.
For non-polar particles in a polar medium, $\alpha_i$ is a negative number.
For more information, see
[J. Israelachvili's book, Chapter 5.](https://www.sciencedirect.com/science/book/9780123751829)

`ionalpha`           | Description
-------------------- | ---------------------------------------
`epsr`               | Relative dielectric constant of medium
`alphaneutral=false` | Do not polarize charged particle

**Limitations:**
Charge-polarizability products for each pair of species is evaluated once during
construction and based on the defined atom types.
Also, `alphaneutral` must be the same for all instances of the potential.
{: .notice--info}


### Cosine Attraction
An attractive potential used for [coarse grained lipids](http://dx.doi.org/10/chqzjk) and with the form,

$$
    \beta u(r) = -\epsilon \cos^2 \left ( \frac{\pi(r-r_c)}{2w_c} \right )
$$

for $r_c\leq r \leq r_c+w_c$. For $r<r_c$, $\beta u=-\epsilon$,
while zero for $r>r_c+w_c$.

`cos2` | Description
------ | --------------------------
`eps`  | Depth, $\epsilon$ (kJ/mol)
`rc`   | Width, $r_c$ (Å)
`wc`   | Decay range, $w_c$ (Å)

### Hard Sphere
`hardsphere`

The hard sphere potential does not take any input. Radii are read from the atomlist at the beginning of the simulation.


### Lennard-Jones

`lennardjones` | Description
-------------- | -------------------------------------------
`mixing=LB`    | Mixing rule. `LB`
`custom`       | Custom $\epsilon$ and $\sigma$ combinations

The Lennard-Jones potential consists of a repulsive and attractive term,

$$
u_{ij}^{\text{LJ}} = 4\epsilon_{ij} \left (
    \left ( \frac{\sigma_{ij}} {r_{ij})} \right )^{12} - \left ( \frac{\sigma_{ij}}{r_{ij})}\right )^6 \right )
$$

where the default mixing rule is Lorentz-Berthelot (`LB`):

$$
\sigma_{ij} = \frac{\sigma_i+\sigma_j}{2} \quad \textrm{and} \quad
\epsilon_{ij} = \sqrt{\epsilon_i \epsilon_j}
$$

The mixing rule can be overridden for specific pairs of atoms:

~~~ yaml
lennardjones:
    custom:
        "Na Cl": {eps: 0.2, sigma: 2}
        "K Cl": {eps: 0.1, sigma: 3}
~~~

### Weeks-Chandler-Andersen

Like Lennard-Jones but cut and shifted to zero at the minimum, forming a
purely [repulsive potential](http://dx.doi.org/ct4kh9),

$$
    u_{ij} = u_{ij}^{\text{LJ}} + \epsilon_{ij} \quad \textrm{for} \quad r<2^{1/6}\sigma_{ij}
$$

`wca`        | Description
------------ | -------------------------------------------
`mixing=LB`  | Mixing rule; only `LB` available.
`custom`     | Custom $\epsilon$ and $\sigma$ combinations


## Bonded Interactions

Bonds and angular potentials are added via the keyword `bondlist` either directly
in a molecule definition (topology) or in energy/bonded where the latter can be
used to add inter-molecular bonds:

~~~ yaml
energy:
    - bonded:
        bondlist: # absolute index
           - harmonic: { index: [56,921], k: 10, req: 15 }
moleculelist:
    - water:
        structure: water.xyz
        bondlist: # index relative to molecule
            - harmonic: { index: [0,1], k: 100, req: 1.0 }  
            - harmonic: { index: [0,2], k: 100, req: 1.0 }
~~~

Bonded potential types:

**Note:**
$\mu V T$ ensembles and Widom insertion are currently unsupported for molecules with bonds.
{: .notice--info}

### Harmonic

`harmonic`     | Harmonic bond
-------------- | -------------------------------------------
`k`            | Harmonic spring constant (kJ/mol/Å$^2$)
`req`          | Equilibrium distance (Å)
`index`        | Array with _exactly two_ index (relative to molecule)

$$
u(r) = \frac{1}{2}k(r-r_{eq})^2
$$

### Finite Extensible Nonlinear Elastic

`fene`         | [Finite Extensible Nonlinear Elastic Potential](http://doi.org/c78x6m)
-------------- | ----------------------------------------------------------------------
`k`            | Bond stiffness (kJ/mol/Å$^2$)
`rmax`         | Maximum separation, $r_m$ (Å)
`index`        | Array with _exactly two_ index (relative to molecule)

$$
u(r) = -\frac{1}{2} k r_{max}^2 \ln \left [ 1-(r/r_m)^2 \right ]
$$

for $r < r_m$; infinity otherwise.

## Geometrical Confinement

`confine`    | Confine molecules to a sub-region
------------ | -------------------------------------------
`type`       | Confinement geometry: `sphere`, `cylinder`, or `cuboid`
`molecules`  | List of molecules to confine (names)
`com=false`  | Apply to molecular mass center
`k`          | Harmonic spring constant in kJ/mol or `inf` for infinity

Confines `molecules` in a given region of the simulation container by applying a harmonic potential on
exterior atom positions, $\mathbf{r}_i$:

$$
U = \frac{1}{2} k \sum^{\text{exterior}} f_i
$$

where $f_i$ is a function that depends on the confinement `type`,
and $k$ is a spring constant. The latter
may be _infinite_ which renders the exterior region strictly inaccessible and may evaluate faster than for finite values.
During equilibration it is advised to use a _finite_ spring constant to drive exterior particles inside the region.

**Note:**
Should you insist on equilibrating with $k=\infty$, ensure that displacement parameters are large enough to transport molecules inside the allowed region, or all moves may be rejected. Further, some analysis routines have undefined behavior for configurations with infinite energies.
{: .notice--danger}

Available values for `type` and their additional keywords:

`sphere`        | Confine in sphere
--------------- | -----------------------------
`radius`        | Radius ($a$)
`origo=[0,0,0]` | Center of sphere ($\mathbf{O}$)
`scale=false`   | Scale radius with volume change, $a^{\prime} = a\sqrt[3]{V^{\prime}/V}$
$f_i$           | $\vert\mathbf{r}_i-\mathbf{O}\vert^2-a^2$

The `scale` option will ensure that the confining radius is scaled whenever the simulation
volume is scaled. This could for example be during a virtual volume move (analysis) or
a volume move in the $NPT$ ensemble.

`cylinder`      | Confine in cylinder along $z$-axis
--------------- | ----------------------------------------
`radius`        | Radius ($a$)
`origo=[0,0,*]` | Center of cylinder ($\mathbf{O}$, $z$-value ignored)
$f_i$           | $\vert (\mathbf{r}_i-\mathbf{O})\circ \mathbf{d}\vert^2 - a^2$

where $\mathbf{d}=(1,1,0)$ and $\circ$ is the entrywise (Hadamard) product.

`cuboid`        | Confine in cuboid
--------------- | --------------------------
`low`           | Lower corner $[x,y,z]$
`high`          | Higher corner $[x,y,z]$
$f_i$           | $\sum_{\alpha\in \{x,y,z\} } (\delta r_{i,\alpha})^2$

where $\delta r$ are distances to the confining, cuboidal faces.
Note that the elements of `low` must be smaller than or equal to the corresponding
elements of `high`.

## Solvent Accessible Surface Area

`sasa`       | SASA Transfer Free Energy
------------ | --------------------------------------------
`radius=1.4` | Probe radius for SASA calculation (angstrom)
`molarity=0` | Molar concentration of co-solute

Calculates the free energy contribution due to

1. atomic surface tension
2. co-solute concentration (typically electrolytes)

via a [SASA calculation](http://doi.org/dws4qm) for each atom.
The energy term is:

$$
    U = \sum_i^N A_{\text{sasa},i} \left ( \gamma_i + c_s \varepsilon_{\text{tfe},i} \right )
$$

where $c_s$ is the molar concentration of the co-solute;
$\gamma_i$ is the atomic surface tension; and $\varepsilon_{\text{tfe},i}$ the atomic transfer free energy,
both specified in the atom topology with `tension` and `tfe`, respectively.

## Penalty Function

This is a version of the flat histogram or Wang-Landau sampling method where
an automatically generated bias or penalty function, $f(\mathcal{X}^d)$,
is applied to the system along a one dimensional ($d=1$)
or two dimensional ($d=2$) reaction coordinate, $\mathcal{X}^d$, so that the configurational integral reads,

$$
    Z(\mathcal{X}^d) = e^{-\beta f(\mathcal{X}^d)} \int e^{-\beta \mathcal{H}(\mathcal{R}, \mathcal{X}^d)} d \mathcal{R}.
$$

where $\mathcal{R}$ denotes configurational space at a given $\mathcal{X}$.
For every visit to a state along the coordinate, a small penalty energy, $f_0$, is
added to $f(\mathcal{X}^d)$ until $Z$ is equal for all $\mathcal{X}$.
Thus, during simulation the free energy landscape is flattened, while the
true free energy is simply the negative of the generated bias function,

$$
    \beta A(\mathcal{X}^d) = -\beta f(\mathcal{X}^d) = -\ln\int e^{-\beta \mathcal{H}(\mathcal{R}, \mathcal{X}^d)}  d \mathcal{R}.
$$

To reduce fluctuations, $f_0$ can be periodically reduced (`update`, `scale`) as $f$ converges.
At the end of simulation, the penalty function is saved to disk as an array ($d=1$) or matrix ($d=2$).
Should the penalty function file be available when starting a new simulation, it is automatically loaded
and used as an initial guess.
This can also be used to run simulations with a _constant bias_ by setting $f_0=0$.

Example setup where the $x$ and $y$ positions of atom 0 are penalized to achieve uniform sampling:

~~~ yaml
energy:
- penalty:
    f0: 0.5
    scale: 0.9
    update: 1000
    file: penalty.dat
    coords:
    - atom: {index: 0, property: "x", range: [-2.0,2.0], resolution: 0.1}
    - atom: {index: 0, property: "y", range: [-2.0,2.0], resolution: 0.1}
~~~

Options:

`penalty`     |  Description
------------- | --------------------
`f0`          |  Penalty energy increment ($kT$)
`update`      |  Interval between scaling of `f0`
`scale=0.8`   |  Scaling factor for `f0`
`nodrift=true`|  Suppress energy drift
`quiet=false` |  Set to true to get verbose output
`file`        |  Name of saved/loaded penalty function
`histogram`   |  Name of saved histogram (not required)
`coords`      |  Array of _one or two_ coordinates

The coordinate, $\mathcal{X}$, can be freely composed by one or two
of the following types (via `coord`):

`atom`        | Single atom properties
------------- | ----------------------------------
`index`       | Atom index
`property`    | `x`, `y`, `z`, `charge`
`range`       | Array w. [min:max] value
`resolution`  | Resolution along coordinate

`cmcm`        | Mass-center separation
------------- | -----------------------------------
`index`       | Array w. exactly two molecule index
`range`       | Array w. [min:max] separation
`resolution`  | Resolution along coordinate
`dir=[1,1,1]` | Directions for distance calc.

`system`      | System property
------------- | -----------------------------------
`property`    | 
`volume`      | System volume

### Multiple Walkers with MPI

If compiled with MPI, the master process collects the bias function from all nodes
upon penalty function `update`.
The _average_ is then re-distributed, offering [linear parallellizing](http://dx.doi.org/10/b5pc4m)
of the free energy sampling. It is crucial that the walk in coordinate space differs on the different
nodes, i.e. by specifying a different random number seed; start configuration; or displacement parameter.
File output and input are prefixed with `mpi{rank}.`

**Information:**
Flat histogram methods are commonly attributed to [Wang and Landau (2001)](http://dx.doi.org/10/bbdg7j)
but the idea appears in earlier works, for example by
[Hunter and Reinhardt (1995)](http://dx.doi.org/10/fns6dq) and
[Engkvist and Karlström (1996)](http://dx.doi.org/10/djjk8z).

