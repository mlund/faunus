<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
});
</script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

# Energy

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
    - confine: {type: sphere, radius: 10, molecules: [water]}
    - nonbonded:
        default: # applied to all atoms
        - lennardjones: {mixing: LB}
        - coulomb: {type: plain, epsr: 1, cutoff: 12}
        Na CH:   # overwrite specific atom pairs
        - wca: { mixing: LB }

    - maxenergy: 100
    - ...
~~~

The keyword `maxenergy` can be used to skip further energy evaluation if a term returns a large
energy change (in kT), which will likely lead to rejection.
The default value is _infinity_.

**Note:**
_Energies_ in MC may contain implicit degrees of freedom, _i.e._ be temperature-dependent,
effective potentials. This is inconsequential for sampling
density of states, but care should be taken when interpreting derived functions such as
energies, entropies, pressure etc.

## Infinite and NaN Energies

In case one or more potential energy terms of the system Hamiltonian returns infinite or NaN energies,
a set of conditions exists to evaluate the acceptance of the proposed move:

- always reject if new energy is NaN (i.e. division by zero)
- always accept if energy change is from NaN to finite energy
- always accept if the energy _difference_ is NaN (i.e. from infinity to minus infinity)

**Note:**
These conditions should be carefully considered if equilibrating a system far from equilibrium.


## External Pressure

This adds the following pressure term (see i.e. [Frenkel and Smith, Chapter 5.4](http://doi.org/c7zg))
to the Hamiltonian, appropriate for MC moves in $\ln V$:

$$
    U = PV - k_BT\left ( N + 1 \right ) \ln V
$$

where $N$ is the total number of molecules and atomic species.

`isobaric`   | Description
-------------|-------------------------------------------------------
`P/unit`     | External pressure where unit can be `mM`, `atm`, or `Pa`.


## Nonbonded Interactions

This term loops over pairs of atoms, $i$, and $j$, summing a given pair-wise additive potential, $u_{ij}$,

$$ U = \sum_{i=0}^{N-1}\sum_{j=i+1}^N u_{ij}(\textbf{r}_j-\textbf{r}_i)$$

Using `nonbonded`, potentials can be arbitrarily mixed and customized for specific particle
combinations. `nonbonded_splined` internally _splines_ the combined potential in an interval [`rmin`,`rmax`] determined
by the following policies:

- `rmin` is decreased towards zero until the potential reaches `u_at_rmin=20` kT
- `rmax` is increased until the potential reaches `u_at_rmax=1e-6` kT

If outside the interval, infinity or zero is returned, respectively.
Finally, the spline precision can be controlled with `utol=1e-5` kT.

Below is a description of possible nonbonded methods. For simple potentials, the hard coded
variants are often the fastest option. For better performance, it is recommended to use `nonbonded_splined` in place of the more robust `nonbonded` method. To check that the combined potential is splined correctly, set `to_disk=true` to print to `A-B_tabulated.dat` the exact and splined combined potentials between species A and B.

`energy`               | $u_{ij}$
---------------------- | ------------------------------------------------------
`nonbonded`            | Any combination of pair potentials (slower, but exact)
`nonbonded_exact`      | An alias for `nonbonded`
`nonbonded_splined`    | Any combination of pair potentials (splined)
`nonbonded_cached`     | Any combination of pair potentials (splined, only intergroup!)
`nonbonded_coulomblj`  | `coulomb`+`lennardjones` (hard coded)
`nonbonded_coulombwca` | `coulomb`+`wca` (hard coded)
`nonbonded_pm`         | `coulomb`+`hardsphere` (fixed `type=plain`, `cutoff`$=\infty$)
`nonbonded_pmwca`      | `coulomb`+`wca` (fixed `type=plain`, `cutoff`$=\infty$)

### Mass Center Cut-offs

For cut-off based pair-potentials working between large molecules, it can be efficient to
use mass center cut-offs between molecular groups, thus skipping all pair-interactions.
A single cut-off can be used between all molecules (`default`), or specified for specific
combinations:

~~~ yaml
- nonbonded:
    cutoff_g2g:
      default: 40
      protein water: 60
~~~

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

coulomb types                            | Keywords      | $\mathcal{S}(q)$
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
U_{\text{reciprocal}} = \frac{2\pi f}{V} \sum_{ {\bf k} \ne {\bf 0}} A_k \vert Q^{q\mu} \vert^2 
$$

$$
U_{\text{surface}} = \frac{ 2\pi f }{ (2\varepsilon_{surf} + 1) V }
\left(
|\sum_{j}q_j{\bf r}_j|^2 + 2 \sum_j q_i {\bf r}_j \cdot \sum_j \boldsymbol{\mu}_j + \vert \sum_j \boldsymbol{\mu}_j \vert^2
\right )
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
{\bf k} = 2\pi\left( \frac{n_x}{L_x} , \frac{n_y}{L_y} ,\frac{n_z}{L_z} \right), {\bf n} \in \mathbb{Z}^3
$$

Like many other electrostatic methods, the Ewald scheme also adds a self-energy term, please see separate Table.
In the case of isotropic periodic boundaries (`ipbc=true`), the orientational degeneracy of the
periodic unit cell is exploited to mimic an isotropic environment, reducing the number
of wave-vectors by one fourth compared with PBC Ewald.
For point charges, [IPBC](http://doi.org/css8) introduce the modification,

$$
Q^q = \sum_j q_j \prod\_{\alpha\in\{x,y,z\}} \cos \left( \frac{2\pi}{L\_{\alpha}} n\_{\alpha} r\_{\alpha,j} \right)
$$

while for point dipoles (currently unavailable),

$$
Q^{\mu} = \sum\_j \boldsymbol{\mu}\_j
\cdot \nabla\_j
\left( \prod\_{ \alpha \in \{ x,y,z \} } \cos \left ( \frac{2\pi}{L\_{\alpha}} n\_{\alpha} r\_{\alpha,j} \right ) \right )
$$

**Limitations:** Ewald summation requires a constant number of particles, i.e. $\mu V T$ ensembles
and Widom insertion are currently unsupported.


### Mean-Field Correction

For cuboidal slit geometries, a correcting mean-field, [external potential](http://dx.doi.org/10/dhb9mj),
$\varphi(z)$, from charges outside the box can be iteratively generated by averaging the charge density, $\rho(z)$, in
$dz$-thick slices along $z$.
This correction assumes that all charges interact with a plain Coulomb potential and that a cubic cutoff is used
via the minimum image convention.

To enable the correction, use the `akesson` keyword at the top level of `energy`:

`akesson`         | Keywords
----------------- | ------------------------------------------------------------
`molecules`       | Array of molecules to operate on
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
used. Through the C++ API or the custom potential explained below, it is easy to add new potentials.

### Charge-Nonpolar

The energy when the field from a point charge, $z_i$, induces a dipole in a polarizable particle of unit-less excess polarizability, $\alpha_j=\left ( \frac{\epsilon_j-\epsilon_r}{\epsilon_j+2\epsilon_r}\right ) a_j^3$, is

$$
    \beta u_{ij} = -\frac{\lambda_B z_i^2 \alpha_j}{2r_{ij}^4}
$$

where $a_j$ is the radius of the non-polar particle and $\alpha_j$ is set in
the atom topology, `alphax`.
For non-polar particles in a polar medium, $\alpha_i$ is a negative number.
For more information, see
[J. Israelachvili's book, Chapter 5.](https://www.sciencedirect.com/science/book/9780123751829)

`ionalpha`           | Description
-------------------- | ---------------------------------------
`epsr`               | Relative dielectric constant of medium

**Limitations:**
Charge-polarizability products for each pair of species is evaluated once during
construction and based on the defined atom types.

### Cosine Attraction

An attractive potential used for [coarse grained lipids](http://dx.doi.org/10/chqzjk) and with the form,

$$
    \beta u(r) = -\epsilon \cos^2 \left ( \frac{\pi(r-r_c)}{2w_c} \right )
$$

for $r\_c\leq r \leq r\_c+w\_c$. For $r<r\_c$, $\beta u=-\epsilon$,
while zero for $r>r\_c+w\_c$.

`cos2` | Description
------ | --------------------------
`eps`  | Depth, $\epsilon$ (kJ/mol)
`rc`   | Width, $r_c$ (Å)
`wc`   | Decay range, $w_c$ (Å)

### Assorted Short Ranged Potentials

The potentials below are often used to keep particles apart and/or to introduce stickiness.
The atomic interaction parameters, e.g., $\sigma_i$ and $\epsilon\_i$, are taken from the
topology.

type             | atomic parameters | $u(r)$ (non-zero part)
---------------- | ----------------- | --------------------------------------------------------
`hardsphere`     | `sigma`           | $\infty$ for $r < \sigma\_{ij}$
`hertz`          | `sigma`, `eps`    | $\epsilon\_{ij} \left ( 1-r / \sigma\_{ij}\right )^{5/2}$ for $r<\sigma\_{ij}$
`lennardjones`   | `sigma`, `eps`    | $4\epsilon\_{ij} \left ( (\sigma\_{ij}/r\_{ij})^{12} - (\sigma\_{ij}/r\_{ij})^6\right )$
`squarewell`     | `sigma`, `eps`    | $-\epsilon\_{ij}$ for $r<\sigma\_{ij}$
[`wca`](http://dx.doi.org/ct4kh9) | `sigma`, `eps` | $u\_{ij}^{\text{LJ}} + \epsilon\_{ij}$ for $r < 2^{1/6}\sigma\_{ij}$

If several potentials are used together and different values for the coefficients are desired,
an aliasing of the parameters' names can be introduced. For example by specifying `sigma: sigma_hs`,
the potential uses the atomic value `sigma_hs` instead of `sigma`, as shown in example below.
To avoid possible conflicts of parameters' names with future keywords of Faunus, we recommend
following naming scheme: `property_pot`, where `property` is either `sigma` or `eps` and
`pot` stands for the potential abbreviation, i.e, `hs`, `hz`, `lj`, `sw`, and `wca`.

Mixing (combination) rules can be specified to automatically parametrize heterogeneous interactions.
If not described otherwise, the same rule is applied to all atomic parameters used by the potential.
No meaningful defaults are defined yet, hence always specify the mixing rule explicitly, e.g.,
`arithmetic` for `hardsphere`.

rule                  | description        | formula
--------------------- | ------------------ | ------------------------------------------------------
`arithmetic`          | arithmetic mean    | $a\_{ij} = \frac 12 \left( a\_{ii} + a_{jj} \right)$
`geometric`           | geometric mean     | $a\_{ij} = \sqrt{a\_{ii} a_{jj}}$
`lorentz_berthelot`   | Lorentz-Berthelot  | `arithmetic` for `sigma`, `geometric` for `eps`

For convenience, the abbreviation `LB` can be used instead of `lorentz_berthelot`.

Custom parameter values can be specified to override the mixing rule for a given pair,
as shown in the example bellow.

~~~ yaml
- lennardjones:
    mixing: LB
    custom:
      Na Cl: {eps: 0.2, sigma: 2}
      K Cl: { ... }

- hertz:
    mixing: LB
    eps: eps_hz
    custom:
      Na Cl: {eps_hz: 0.2, sigma: 2}

- hardsphere:
    mixing: arithmetic
    sigma: sigma_hs
    custom:
      Na Cl: {sigma_hs: 2}
~~~

### SASA (pair potential)

This calculates the surface area of two intersecting particles or radii $R$ and $r$
to estimate an energy based on transfer-free-energies (TFE) and surface tension.
The total surface area is calculated as

$$
    A = 4\pi \left ( R^2 + r^2 \right ) - 2\pi \left (  Rh_1 + rh_2  \right )
$$

where $h\_1$ and $h\_2$ are the heights of the spherical caps comprising the lens
formed by the overlapping spheres. For complete overlap, or when far apart, the
full area of the bigger sphere or the sum of both spheres are returned.
The pair-energy is calculated as:

$$
    u_{ij} = A \left ( \gamma_{ij} + c_s \varepsilon_{\text{tfe},ij} \right )
$$

where $\gamma\_{ij}$ and $\varepsilon\_{\text{tfe},ij}$ are the arithmetic means of
 `tension` and `tfe` provided in the atomlist.

Note that SASA is strictly not additive and this pair-potential is merely
a poor-mans way of approximately take into account ion-specificity and
hydrophobic/hydrophilic interactions. Faunus offers also a full, albeit yet
experimental implementation of [Solvent Accessible Surface Area] energy.

`sasa`       | Description
------------ | ----------------------------------------------------------
`molarity`   | Molar concentration of co-solute, $c_s$
`radius=1.4` | Probe radius for SASA calculation (angstrom)
`shift=true` | Shift to zero at large separations

### Custom

This takes a user-defined expression and a list of constants to produce a runtime,
custom pair-potential.
While perhaps not as computationally efficient as hard-coded potentials, it is a
convenient way to access alien potentials. Used in combination with `nonbonded_splined`
there is no overhead since all potentials are splined.

`custom`     | Description
------------ | --------------------------------------------------------
`function`   | Mathematical expression for the potential (units of kT)
`constants`  | User-defined constants
`cutoff`     | Spherical cut-off distance

The following illustrates how to define a Yukawa potential:

~~~ yaml
custom:
    function: lB * q1 * q2 / r * exp( -r/D ) # in kT
    constants:
        lB: 7.1  # Bjerrum length
        D: 30    # Debye length
~~~

The function is passed using the efficient
[ExprTk library](http://www.partow.net/programming/exprtk/index.html) and
a rich set of mathematical functions and logic is available.
In addition to user-defined constants, the following symbols are defined:

`symbol`   | Description
---------- | ---------------------------------------------
`e0`       | Vacuum permittivity [C^2/J/m]
`inf`      | infinity
`kB`       | Boltzmann's constant [J/K]
`kT`       | Boltzmann's constant x temperature [J]
`Nav`      | Avogadro's number [1/mol]
`pi`       | Pi
`q1`,`q2`  | particle charges [e]
`r`        | particle-particle separation [angstrom]
`Rc`       | Spherical cut-off [angstrom]
`s1`,`s2`  | particle sigma [angstrom]
`T`        | temperature [K]

## Custom External Potential

This applies a custom external potential to atoms or molecular mass centra
using the [ExprTk library](http://www.partow.net/programming/exprtk/index.html)
syntax.

`customexternal` | Description
---------------- | --------------------------------------------------------
`molecules`      | Array of molecules to operate on
`com=false`      | Operate on mass-center instead of individual atoms?
`function`       | Mathematical expression for the potential (units of kT)
`constants`      | User-defined constants

In addition to user-defined `constants`, the following symbols are available:

`symbol`   | Description
---------- | ---------------------------------------
`e0`       | Vacuum permittivity [C^2/J/m]
`inf`      | infinity
`kB`       | Boltzmann's constant [J/K]
`kT`       | Boltzmann's constant x temperature [J]
`Nav`      | Avogadro's number [1/mol]
`pi`       | Pi
`q`        | particle charge [e]
`s`        | particle sigma [angstrom]
`x`,`y`,`z`| particle positions [angstrom]
`T`        | temperature [K]

If `com=true`, charge refers to the molecular net-charge, and `x,y,z` the mass-center coordinates.
The following illustrates how to confine molecules in a spherical shell of radius, _r_, and
thickness _dr_:

~~~ yaml
customexternal:
    molecules: [water]
    com: true
    constants: {radius: 15, dr: 3}
    function: >
        var r2 := x^2 + y^2 + z^2;
        if ( r2 < radius^2 )
           1000 * ( radius-sqrt(r2) )^2;
        else if ( r2 > (radius+dr)^2 )
           1000 * ( radius+dr-sqrt(r2) )^2;
        else
           0;
~~~


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
        structure: "water.xyz"
        bondlist: # index relative to molecule
            - harmonic: { index: [0,1], k: 5024, req: 0.9572 }
            - harmonic: { index: [0,2], k: 5024, req: 0.9572 }
            - harmonic_torsion: { index: [1,0,2], k: 628, aeq: 104.52 }
~~~

Bonded potential types:

**Note:**
$\mu V T$ ensembles and Widom insertion are currently unsupported for molecules with bonds.

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
`index`        | Array with _exactly two_ indices (relative to molecule)

Finite extensible nonlinear elastic potential long range repulsive potential.

$$
     u(r) =
  \begin{cases} 
   -\frac{1}{2} k r_{\mathrm{max}}^2 \ln \left [ 1-(r/r_{\mathrm{max}})^2 \right ],       & \text{if } r < r_{\mathrm{max}} \\
   \infty, & \text{if } r \geq r_{\mathrm{max}}
  \end{cases}
$$

**Note:**
It is recommend to only use the potential if the initial configuration is near equilibrium, which prevalently depends on the value of `rmax`.
Should one insist on conducting simulations far from equilibrium, a large displacement parameter is recommended to reach finite energies.

### Finite Extensible Nonlinear Elastic + WCA

`fene+wca`     | [Finite Extensible Nonlinear Elastic Potential + WCA](http://doi.org/c78x6m)
-------------- | ----------------------------------------------------------------------
`k`            | Bond stiffness (kJ/mol/Å²)
`rmax`         | Maximum separation, $r_m$ (Å)
`eps=0`        | Epsilon energy scaling (kJ/mol)
`sigma=0`      | Particle diameter (Å)
`index`        | Array with _exactly two_ indices (relative to molecule)

Finite extensible nonlinear elastic potential long range repulsive potential combined
with the short ranged Weeks-Chandler-Anderson (wca) repulsive potential. This potential is particularly useful in combination with the `nonbonded_cached` non-bonded energy.

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

Note that the implementation of Solvent Accessible Surface Area potential is considered _experimental_.
The code is untested, unoptimized and the configuration syntax below can change.
The FreeSASA library option has to be enabled when [compiling].

`sasa`       | SASA Transfer Free Energy
------------ | --------------------------------------------
`radius=1.4` | Probe radius for SASA calculation (angstrom)
`molarity`   | Molar concentration of co-solute

Calculates the free energy contribution due to

1. atomic surface tension
2. co-solute concentration (typically electrolytes)

via a [SASA calculation](http://dx.doi.org/10/dbjh) for each atom, as implemented in
the [FreeSASA library](https://freesasa.github.io/).

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

Flat histogram methods are often attributed to [Wang and Landau (2001)](http://dx.doi.org/10/bbdg7j)
but the idea appears in earlier works, for example by
[Hunter and Reinhardt (1995)](http://dx.doi.org/10/fns6dq) and
[Engkvist and Karlström (1996)](http://dx.doi.org/10/djjk8z).

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

`penalty`        |  Description
---------------- | --------------------
`f0`             |  Penalty energy increment (_kT_)
`update`         |  Interval between scaling of `f0`
`scale`          |  Scaling factor for `f0`
`nodrift=true`   |  Suppress energy drift
`quiet=false`    |  Set to true to get verbose output
`file`           |  Name of saved/loaded penalty function
`overwrite=true` |  If `false`, don't save final penalty function
`histogram`      |  Name of saved histogram (not required)
`coords`         |  Array of _one or two_ coordinates

The coordinate, $\mathcal{X}$, can be freely composed by one or two
of the types listed in the next section (via `coords`).


### Reaction Coordinates

The following reaction coordinates can be used for penalising the energy and can further
be used when analysing the system (see Analysis).

#### Atom Properties

`coords=[atom]`| Single atom properties
-------------- | ----------------------------------
`index`        | Atom index
`property`     | `x`, `y`, `z`, `q`
`range`        | Array w. [min:max] value
`resolution`   | Resolution along coordinate
`N`            | Number of atoms of ID=`index`

#### Molecule Properties

`coords=[molecule]`       | Single molecule properties
------------------------- | --------------------------------------------------------------------------------------------------------------------------------
`index`                   | Molecule index
`range`                   | Array w. [min:max] value
`resolution`              | Resolution along coordinate
`property`                | Options:
`angle`                   | Angle between instantaneous principal axis and given `dir` vector
`com_x`, `com_y`, `com_z` | Mass center coordinates
`confid`                  | Conformation id corresponding to frame in `traj` (see molecular topology).
`end2end`                 | Distance between first and last atom
`mu_x`, `mu_y`, `mu_z`    | Molecular dipole moment components
`mu`                      | Molecular dipole moment scalar (eA/charge)
`muangle`                 | Angle between dipole moment and given `dir` vector
`N`                       | Number of atoms in group
`Q`                       | Monopole moment (net charge)
`atomatom`                | Distance along `dir` between 2 atoms specified by the `indexes` array
`cmcm`                    | Absolute mass-center separation between groups defined by the intervals `indexes[0]:indexes[1]` and `indexes[2]:indexes[3]`
`cmcm_z`                  | _z_-component of mass-center separation between groups defined by the intervals `indexes[0]:indexes[1]` and `indexes[2]:indexes[3]`
`L/R`                     | Ratio between height and radius of a cylindrical lipid vesicle (ad-hoc RC for bending modulus calculations)

Notes:

- the molecular dipole moment is defined w. respect to the mass-center
- for `angle`, the principal axis is the eigenvector corresponding to the smallest eigenvalue of the gyration tensor

#### System Properties

`coords=[system]` | System property
----------------- | ----------------------------------------
`range`           | Array w. [min:max] value
`resolution`      | Resolution along coordinate
`property`        | Options:
`V`               | System volume
`Q`               | System net-charge
`Lx`,`Ly`,`Lz`    | Side lengths of enclosing cuboid
`height`          | Alias for `Lz`
`radius`          | Radius of spherical or cylindrical geometries
`N`               | Number of active particles

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

The following starts all MPI processes with the same input file and MPI prefix is automatically
appended to all other input and output:

~~~ bash
yason.py input.yml | mpirun --np 6 --stdin all faunus -s state.json
~~~

Here, each process automatically looks for `mpi{nproc}.state.json`.

## Constraining the system

Reaction coordinates can be used to constrain the system within a `range`
using the `constrain` energy term. Stepping outside the range results in an inifinite
energy, forcing rejection. For example,

~~~ yaml
energy:
    - constrain: {type: molecule, index: 0, property: end2end, range: [0,200]}
~~~

Tip: placing `constrain` at the _top_ of the energy list is more efficient as the remaining
energy terms are skipped should an infinite energy arise.

