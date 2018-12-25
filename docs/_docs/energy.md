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

The energy terms are specified in `energy` at the top level input and evaluated in the order given.
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

    - maxenergy: 100
    - ...
~~~

The keyword `maxenergy` can be used to stop energy evaluation when a large
energy change (in kT) is observed which will most likely lead to rejection.
The default value is _infinity_.

**Note:**
_Energies_ in MC may contain implicit degrees of freedom, _i.e._ be temperature-dependent,
effective potentials. This is inconsequential for sampling
density of states, but care should be taken when interpreting derived functions such as
energies, entropies, pressure etc.
{: .notice--info}

## Infinite and NaN Energies

In case one or more potential energy terms of the system Hamiltonian returns infinite or NaN energies,
a set of conditions exists to evaluate the acceptance of the proposed move:

- always reject if new energy is NaN (i.e. division by zero)
- always accept if energy change is from NaN to finite energy
- always accept if the energy _difference_ is NaN (i.e. from infinity to minus infinity)

**Note:**
These conditions should be carefully considered if equilibrating a system far from equilibrium,
particularly if using discontinuous potentials.
{: .notice--notice}


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
`nonbonded_coulombwca` | `coulomb`+`wca`
`nonbonded_pm`         | `coulomb`+`hardsphere` (fixed `type=plain`, `cutoff`$=\infty$)
`nonbonded_pmwca`      | `coulomb`+`wca` (fixed `type=plain`, `cutoff`$=\infty$)


### OpenMP Control

If compiled with OpenMP, the following keywords can be used to control parallelisation
for non-bonded interactions. The best combination depends on the simulated system size and
composition. Currently, parallelisation is disabled by default.

~~~ yaml
- nonbonded:
    openmp: [g2g, i2all]
~~~

`openmp`  | Description
--------- | -------------------------------------------
`g2g`     | Distribute on a molecule-to-molecule basis 
`i2all`   | Parallelise single particle energy evaluations


## Electrostatics

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

### Ewald Summation

If type is `ewald`, terms from reciprocal space; surface energies; and
self energies are automatically added to the Hamiltonian, activating additional keywords:

`type=ewald`         | Description
-------------------- | ---------------------------------------------------------------------
`kcutoff`            | Reciprocal-space cutoff
`epss=0`             | Dielectric constant of surroundings, $\varepsilon_{surf}$ (0=tinfoil)
`ipbc=false`         | Use isotropic periodic boundary conditions, [IPBC](http://doi.org/css8).
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
For point charges, [IPBC](http://doi.org/css8) introduce the modification,

$$
Q^q = \sum_jq_j\prod_{\alpha\in\{x,y,z\}}\cos\left(\frac{2\pi}{L_{\alpha}}n_{\alpha} r_{\alpha,j}\right)
$$

while for point dipoles (currently unimplemented),

$$
Q^{\mu} = \sum_j\boldsymbol{\mu}_j\cdot\nabla_j\left(\prod_{\alpha \in\{x,y,z\}}\cos\left(\frac{2\pi}{L_{\alpha}}n_{\alpha}r_{\alpha,j}\right)\right).
$$

**Limitations:** Ewald summation requires a constant number of particles, i.e. $\mu V T$ ensembles
and Widom insertion are currently unsupported.
{: .notice--info}

### Mean-Field Correction

For cuboidal slit geometries, a correcting mean-field, [external potential](http://dx.doi.org/10/dhb9mj),
$\varphi(z)$, from charges outside the box can be iteratively generated by averaging the charge density, $\rho(z)$, in
$dz$-thick slices along $z$.
This correction assumes that all charges interact with a plain Coulomb potential and that a cubic cutoff is used
via the minimum image convention.

To enable the correction, use the `akesson` keyword at the top level of `energy`:

`akesson`         | Keywords
----------------- | ------------------------------------------------------------
`epsr`            | Relative dielectric constant
`nstep`           | Number of energy evalutations between updating $\rho(z)$
`dz=0.2`          | $z$ resolution (angstrom)
`nphi=10`         | Multiple of `nstep` in between updating $\varphi(z)$
`file=mfcorr.dat` | File with $\rho(z)$ to either load or save
`fixed=false`     | If true, assume that `file` is converged. No further updating and faster.

The density is updated every `nstep` energy calls, while the external potential can be updated
slower (`nphi`) since it affects the ensemble.
A reasonable value of `nstep` is system dependent and can be a rather large value.
Updating the external potential on the fly leads to energy drifts that decrease for consecutive runs.
Production runs should always be performed with `fixed=true` and a well converged $\rho(z)$.

At the end of simulation, `file` is overwritten unless `fixed=true`.

## Pair Potentials

In addition to the Coulombic pair-potentials described above, a number of other pair-potentials can be
used. Through the C++ API, it is easy to add new potentials.

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

### SASA Pair-potential

This calculates the surface area of two intersecting particles or radii $R$ and $r$
to estimate an energy based on transfer-free-energies (TFE) and surface tension.
The total surface area is calculated as

$$
    A = 4\pi \left ( R^2 + r^2 \right ) - 2\pi \left (  Rh_1 + rh_2  \right )
$$

where $h_1$ and $h_2$ are the heights of the spherical caps comprising the lens
formed by the overlapping spheres. For complete overlap, or when far apart, the
full area of the bigger sphere or the sum of both spheres are returned.
The pair-energy is calculated as:

$$
    u_{ij} = A \left ( \gamma_{ij} + c_s \varepsilon_{\text{tfe},ij} \right )
$$

where $\gamma_{ij}$ and $\varepsilon_{\text{tfe},ij}$ are the arithmetic means of
 `tension` and `tfe` provided in the atomlist.

Note that SASA is strictly not additive and this
pair-potential is merely a poor-mans way of approximately take into account ion-specificity
and hydrophobic/hydrophilic interactions.

`sasa`       | Description
------------ | ----------------------------------------------------------
`molarity`   | Molar concentration of co-solute, $c_s$.
`radius=1.4` | Probe radius for SASA calculation (angstrom)
`shift=true` | Shift to zero at large separations


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
    - water: # TIP3P
        structure: water.xyz
        bondlist: # index relative to molecule
            - harmonic: { index: [0,1], k: 5024, req: 0.9572 }
            - harmonic: { index: [0,2], k: 5024, req: 0.9572 }
            - harmonic_torsion: { index: [1,0,2], k: 628, aeq: 104.52 }
~~~

Bonded potential types:

**Note:**
$\mu V T$ ensembles and Widom insertion are currently unsupported for molecules with bonds.
{: .notice--info}

### Harmonic

`harmonic`     | Harmonic bond
-------------- | -------------------------------------------
`k`            | Harmonic spring constant (kJ/mol/Å²)
`req`          | Equilibrium distance (Å)
`index`        | Array with _exactly two_ indices (relative to molecule)

$$
u(r) = \frac{1}{2}k(r-r_\mathrm{eq})^2
$$

### Finite Extensible Nonlinear Elastic

`fene`         | [Finite Extensible Nonlinear Elastic Potential](http://doi.org/c78x6m)
-------------- | ----------------------------------------------------------------------
`k`            | Bond stiffness (kJ/mol/Å²)
`rmax`         | Maximum separation, $r_m$ (Å)
`eps=0`        | Epsilon energy scaling (kJ/mol)
`sigma=0`      | Particle diameter (Å)
`index`        | Array with _exactly two_ indices (relative to molecule)

Finite extensible nonlinear elastic potential long range repulsive potential combined
with the short ranged Weeks-Chandler-Anderson (wca) repulsive potential.

$$
     u(r) =
  \begin{cases} 
   -\frac{1}{2} k r_{\mathrm{max}}^2 \ln \left [ 1-(r/r_{\mathrm{max}})^2 \right ] + u_{\mathrm{wca}}, & \text{if } 0 < r \leq 2^{1/6}\sigma \\
   -\frac{1}{2} k r_{\mathrm{max}}^2 \ln \left [ 1-(r/r_{\mathrm{max}})^2 \right ],       & \text{if } 2^{1/6}\sigma < r < r_{\mathrm{max}} \\
   \infty, & \text{if } r \geq r_{\mathrm{max}}
  \end{cases}
$$

**Note:**
It is recommend to only use the potential if the initial configuration is near equilibrium, which prevalently depends on the value of `rmax`.
Should one insist on conducting simulations far from equilibrium, a large displacement parameter is recommended to reach finite energies.
{: .notice--info}


### Harmonic torsion

`harmonic_torsion` | Harmonic torsion
------------------ | -------------------------------------------
`k`                | Harmonic spring constant (kJ/mol/rad²)
`aeq`              | Equilibrium angle $\alpha_\mathrm{eq}$ (deg)
`index`            | Array with _exactly three_ indices (relative to molecule)

$$
u(r) = \frac{1}{2}k(\alpha - \alpha_\mathrm{eq})^2
$$

### Cosine based torsion (GROMOS-96)

`g96_torsion`      | Cosine based torsion
------------------ | -------------------------------------------
`k`                | Force constant (kJ/mol)
`aeq`              | Equilibrium angle $\alpha_\mathrm{eq}$ (deg)
`index`            | Array with _exactly three_ indices (relative to molecule)

$$
u(r) = \frac{1}{2}k(\cos(\alpha) - \cos(\alpha_\mathrm{eq}))^2
$$

### Proper periodic dihedral

`periodic_dihedral` | Proper periodic dihedral
------------------- | -------------------------------------------
`k`                 | Force constant (kJ/mol)
`n`                 | Periodicity (multiplicity) of the dihedral (integer)
`phi`               | Angle $\phi_\mathrm{syn}$ (deg)
`index`             | Array with _exactly four_ indices (relative to molecule)

$$
u(r) = k(1 + \cos(n\phi - \phi_\mathrm{syn}))
$$


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
`molarity`   | Molar concentration of co-solute

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

Example setup where the charges of atom 0 and 1, $q_0$ and $q_1$, are penalized within the ranges [-0.5, +0.5] and [-1.0, 0.0], respectively.
The method can be used to estimate for example free energy of binding.
The penalty assigned to the charges is carried out as the charges are moved within their respective ranges.

~~~ yaml
energy:
- penalty:
    f0: 0.5
    scale: 0.9
    update: 1000
    file: penalty.dat
    histogram: penalty.dat
    coords:
    - atom: {index: 0, property: "q", range: [-0.5,0.5], resolution: 0.01}
    - atom: {index: 1, property: "q", range: [-1.0,0.0], resolution: 0.01}
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
of the types listed in the next section (via `coords`)

**Note:**
When using penalty energies there may currently be unresolved issues loading
previously saved states. If so, the program will terminate.
{: .notice--info}


### Reaction Coordinates

The following reaction coordinates can be used for penalising the energy and can further
be used when analysing the system (see Analysis).

#### Atom Properties

`atom`        | Single atom properties
------------- | ----------------------------------
`index`       | Atom index
`property`    | `x`, `y`, `z`, `q`
`range`       | Array w. [min:max] value
`resolution`  | Resolution along coordinate

#### Molecule Properties

`molecule`                | Single molecule properties
------------------------- | ----------------------------------
`index`                   | Molecule index
`range`                   | Array w. [min:max] value
`resolution`              | Resolution along coordinate
`property`                |
`angle`                   | Angle between instantaneous principal axis and given `dir` vector
`com_x`, `com_y`, `com_z` | Mass center coordinates
`confid`                  | Conformation id corresponding to frame in `traj` (see molecular topology).
`mu_x`, `mu_y`, `mu_z`    | Molecular dipole moment components
`mu`                      | Molecular dipole moment scalar (eA/charge)
`muangle`                 | Angle between dipole moment and given `dir` vector
`N`                       | Number of atoms in group
`Q`                       | Monopole moment (net charge)

Notes:

- the molecular dipole moment is defined w. respect to the mass-center
- for `angle`, the principle axis is calculated by diagonalising the gyration tensor

#### Molecule Separation

This returns the minimum distance between the mass centers of two molecules.
Useful for calculating i.e. the potential of mean force between strongly
interacting molecular groups.

`cmcm`        | Mass-center separation
------------- | -----------------------------------
`index`       | Array w. exactly two molecule index
`range`       | Array w. [min:max] separation
`resolution`  | Resolution along coordinate
`dir=[1,1,1]` | Directions for distance calc.

#### System Properties

`system`       | System property
-------------- | ----------------------------------------------
`range`        | Array w. [min:max] value
`resolution`   | Resolution along coordinate
`property`     |
`V`            | System volume
`Q`            | System net-charge
`Lx`,`Ly`,`Lz` | Side lengths of enclosing cuboid
`height`       | Alias for `Lz`
`radius`       | Radius of spherical or cylindrical geometries

The enclosing cuboid is the smallest cuboid that can contain the geometry.
For example, for a cylindrical simulation container, `Lz` is the height
and `Lx=Ly` is the diameter.


### Multiple Walkers with MPI

If compiled with MPI, the master process collects the bias function from all nodes
upon penalty function `update`.
The _average_ is then re-distributed, offering [linear parallellizing](http://dx.doi.org/10/b5pc4m)
of the free energy sampling. It is crucial that the walk in coordinate space differs on the different
nodes, i.e. by specifying a different random number seed; start configuration; or displacement parameter.
File output and input are prefixed with `mpi{rank}.`

**Information:**
Flat histogram methods are often attributed to [Wang and Landau (2001)](http://dx.doi.org/10/bbdg7j)
but the idea appears in earlier works, for example by
[Hunter and Reinhardt (1995)](http://dx.doi.org/10/fns6dq) and
[Engkvist and Karlström (1996)](http://dx.doi.org/10/djjk8z).

