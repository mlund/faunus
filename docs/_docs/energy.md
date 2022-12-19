<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
});
</script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

# Energy

The system energy, or Hamiltonian, consists of a sum of potential energy terms,

$$
\mathcal{H}\_{sys} = U\_1 + U\_2 + ...
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
        - coulomb: {type: plain, epsr: 1}
        Na CH:   # overwrite specific atom pairs
        - wca: { mixing: LB }

    - maxenergy: 100
    - ...
~~~

The keyword `maxenergy` can be used to skip further energy evaluation if a term returns a large
energy change (in kT), which will likely lead to rejection.
The default value is _infinity_.

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

These conditions should be carefully considered if equilibrating a system far from equilibrium.


## External Pressure

This adds the following pressure term (see i.e. [Frenkel and Smith, Chapter 5.4](http://doi.org/c7zg))
to the Hamiltonian, appropriate for MC moves in $\ln V$:

$$
    U = PV - k_BT\left ( N + 1 \right ) \ln V
$$

where $N$ is the total number of molecules and atomic species.

`isobaric`   | Description
-------------|-------------------------------------------------------------------
`P/unit`     | External pressure where unit can be `mM`, `atm`, `Pa`, `bar`, `kT`


## Nonbonded Interactions

This term loops over pairs of atoms, $i$, and $j$, summing a given pair-wise additive potential, $u_{ij}$,

$$ U = \sum_{i=0}^{N-1}\sum_{j=i+1}^N u_{ij}(\textbf{r}\_j-\textbf{r}\_i)$$

The most general method is `nonbonded` where potentials can be arbitrarily mixed and customized
for specific particle combinations.

Example:

~~~ yaml
- nonbonded:
    default: # default pair potential
        - lennardjones: {mixing: LB}
        - coulomb: {type: fanourgakis, epsr: 1.0, cutoff: 12}
    Ow Ca: # custom potential for atom type "Ow" and atom type "Ca"
        - wca: {mixing: LB}
~~~

Below is a description of possible nonbonded methods. For simple potentials, the hard coded
variants are often the fastest option.
For better performance, it is recommended to use `nonbonded_splined` in place of the more robust `nonbonded` method.

`energy`               | $u\_{ij}$
---------------------- | ------------------------------------------------------
`nonbonded`            | Any combination of pair potentials (slower, but exact)
`nonbonded_exact`      | An alias for `nonbonded`
`nonbonded_splined`    | Any combination of pair potentials (splined)
`nonbonded_cached`     | Any combination of pair potentials (splined, only intergroup!)
`nonbonded_coulomblj`  | `coulomb`+`lennardjones` (hard coded)
`nonbonded_coulombwca` | `coulomb`+`wca` (hard coded)
`nonbonded_pm`         | `coulomb`+`hardsphere` (fixed `type=plain`, `cutoff`$=\infty$)
`nonbonded_pmwca`      | `coulomb`+`wca` (fixed `type=plain`, `cutoff`$=\infty$)


### Mass Center Cutoffs

For cutoff based pair-potentials working between large molecules, it can be efficient to
use mass center cutoffs between molecular groups, thus skipping all pair-interactions.
A single cutoff can be used between all molecules (`default`), or specified for specific
combinations:

~~~ yaml
- nonbonded:
    cutoff_g2g:
      default: 40.0
      protein polymer: 20.0
~~~

If `default` is omitted, only the specified pairs are subject to the cutoffs.
Finally, `cutoff_g2g: 40.0` is allowed for a uniform cutoff between all groups.


### Spline Options

The `nonbonded_splined` method internally _splines_ the potential in an automatically determined
interval [`rmin`,`rmax`] determined by the following policies:

- `rmin` is decreased towards zero until the potential reaches `u_at_rmin`.
- `rmax` is increased until the potential reaches `u_at_rmax`.

If above the interval, zero is returned.
If below the interval, the exact energy (or infinity) is returned.
For details about the splines for each pair, use
`to_disk` and/or maximize the verbosity level (`--verbosity`) when
running faunus.

Keyword            | Description
------------------ | ------------------------------------------------
`utol=1e-3`        | Spline precision
`u_at_rmin=20`     | Energy threshold at short separations (_kT_)
`u_at_rmax=1e-6`   | Energy threshold at long separations (_kT_)
`to_disk=False`    | Create datafiles w. exact and splined potentials
`hardsphere=False` | Use hardsphere repulsion below rmin

Note: Anisotropic pair-potentials cannot be splined. This also applies
to non-shifted electrostatic potentials such as `plain` and un-shifted `yukawa`.

## Parallel summation

Depending on how Faunus was compiled, parallel, nonbonded summation schemes may be available.
Activate with:

~~~ yaml
- nonbonded:
    summation_policy: parallel
    ...
~~~

where `parallel` uses C++ internal threading; `openmp` uses OpenMP; and `serial` skip
parallel summation (default). A warning will be issued if the desired scheme is unavailable.
For the `openmp` policy, you may control the number of threads with the environmental variable
`OMP_NUM_THREADS`.
Summation policies other than `serial` may require substantial memory for systems with many particles.


## Electrostatics

 `coulomb`             |  Description
 --------------------- |  -------------------------------------------------
 `type`                |  Coulomb type, see below
 `cutoff`              |  Spherical cutoff, $R\_c$ (Å) after which the potential is zero
 `epsr`                |  Relative dielectric constant of the medium
 `utol=0.005/lB`       |  Error tolerence for splining; default value depends on the Bjerrum length, lB
 `debyelength=`$\infty$|  Debye length (Å) if using `ewald`, `poisson`, `yukawa`

This is a multipurpose potential that handles several electrostatic methods.
Beyond a spherical real-space cutoff, $R\_c$, the potential is zero while if
below,

$$
\tilde{u}^{(zz)}\_{ij}(\bar{r}) = \frac{e^2 z\_i z\_j }{ 4\pi\epsilon\_0\epsilon\_r |\bar{r}| }\mathcal{S}(q)
$$

where $\bar{r} = \bar{r}\_j - \bar{r}\_i$, and tilde indicate that a short-range function $\mathcal{S}(q=|\bar{r}|/R\_c)$ is used to trucate the interactions. The available short-range functions are:

coulomb types                            | Keywords          | $\mathcal{S}(q)$
---------------------------------------- | ----------------- | ---------------------------------------------------
[`plain`](http://doi.org/ctnnsj)         |                   | 1
[`ewald`](http://doi.org/dgpdmc)         | `alpha`           | $\frac{1}{2}\text{erfc}\left(\alpha R\_c q + \frac{\kappa}{2\alpha}\right)\text{exp}\left(2\kappa R\_c q\right) + \frac{1}{2}\text{erfc}\left(\alpha R\_c q - \frac{\kappa}{2\alpha}\right)$
[`reactionfield`](http://doi.org/dbs99w) | `epsrf`           | $1+\frac{\epsilon_{RF}-\epsilon\_r}{2\epsilon_{RF}+\epsilon\_r}q^3-3\frac{\epsilon_{RF}}{2\epsilon_{RF}+\epsilon\_r}q$
[`poisson`](http://doi.org/10/c5fr)      | `C=3`, `D=3`      | $(1-\acute{q})^{D+1}\sum_{c=0}^{C-1}\frac{C-c}{C}{D-1+c\choose c}\acute{q}^c$
[`qpotential`](https://arxiv.org/abs/1904.10335) | `order`   | $\prod_{n=1}^{\text{order}}(1-q^n)$
[`fanourgakis`](http://doi.org/f639q5)   |                   | $1-\frac{7}{4}q+\frac{21}{4}q^5-7q^6+\frac{5}{2}q^7$
[`fennell`](http://doi.org/10/bqgmv2)    | `alpha`           | $\text{erfc}(\alpha R\_cq)-q\text{erfc}(\alpha R\_c)+(q-1)q\left(\text{erfc}(\alpha R\_c)+\frac{2\alpha R\_c}{\sqrt{\pi}}\text{exp}(-\alpha^2R\_c^2)\right)$
[`zerodipole`](http://doi.org/10/fhcfn4) | `alpha`           | $\text{erfc}(\alpha R\_cq)-q\text{erfc}(\alpha R\_c)+\frac{(q^2-1)}{2}q\left(\text{erfc}(\alpha R\_c)+\frac{2\alpha R\_c}{\sqrt{\pi}}\text{exp}(-\alpha^2R\_c^2)\right)$
[`zahn`](http://doi.org/10/cmx5vd)       | `alpha`           | $\text{erfc}(\alpha R\_c q)-(q-1)q\left(\text{erfc}(\alpha R\_c)+\frac{2\alpha R\_c}{\sqrt{\pi}}\text{exp}(-\alpha^2R\_c^2)\right)$
[`wolf`](http://doi.org/cfcxdk)          | `alpha`           | $\text{erfc}(\alpha R\_cq)-\text{erfc}(\alpha R\_c)q$
`yukawa`                                 | `debyelength`, `shift=false` | As `plain` with screening or, if shifted, `poisson` with `C=1` and `D=1`

Internally $\mathcal{S}(q)$ is _splined_ whereby all types evaluate at similar speed.
For the `poisson` potential,

$$
\acute{q} = \frac{1-\exp\left(2\kappa R\_c q\right)}{1-\exp\left(2\kappa R\_c\right)}
$$

which as the inverse Debye length, $\kappa\to 0$ gives $\acute{q}=q$.
The `poisson` scheme can generate a number of other truncated pair-potentials found in the litterature,
depending on `C` and `D`. Thus, for an infinite Debye length, the following holds:

`C` | `D` | Equivalent to
--- | --- | ----------------------
 1  | -1  | Plain Coulomb within cutoff, zero outside
 1  |  0  | [Undamped Wolf](http://doi.org/10.1063/1.478738)
 1  |  1  | [Levitt](http://doi.org/10/fp959p) / [Undamped Fenell](http://doi.org/10/bqgmv2)
 1  |  2  | [Kale](http://doi.org/10/csh8bg)
 1  |  3  | [McCann](http://doi.org/10.1021/ct300961)
 2  |  1  | [Undamped Fukuda](http://doi.org/10.1063/1.3582791)
 2  |  2  | [Markland](http://doi.org/10.1016/j.cplett.2008.09.019)
 3  |  3  | [Stenqvist](http://doi.org/10/c5fr)
 4  |  3  | [Fanourgakis](http://doi.org/10.1063/1.3216520)


### Debye Screening Length

A background screening due to implicit ions can be added by specifying the keyword `debyelength` to the schemes

- `yukawa`
- `ewald`
- `poisson`

The `yukawa` scheme has simple exponential screening and, like `plain`, an _infinite_ cutoff.
If `shift: true` is passed to the yukawa scheme, the potential is shifted to give zero potential and force 
at the now _finite_ `cutoff` distance (simply an alias for `poisson` with C=1 and D=1).
The list below shows alternative ways to specify the background electrolyte, and will automatically deduce
the salt stoichiometry based on valencies:

~~~ yaml
    debyelength: 30.0, epsr: 79.8       # assuming 1:1 salt, e.g. NaCl
    molarity: 0.02                      # 0.02 M 1:1 salt, e.g. NaCl
    molarity: 0.01, valencies: [2,3,-2] # 0.01 M Ca₂Al₂(SO₄)₅
~~~


### Multipoles

If the type `coulomb` is replaced with `multipole` then the electrostatic energy will in addition to
monopole-monopole interactions include contributions from monopole-dipole, and dipole-dipole
interactions. Multipolar properties of each particle is specified in the Topology.
The `zahn` and `fennell` approaches have undefined dipolar self-energies and are therefore not recommended for such systems.

The ion-dipole interaction is described by

$$
\tilde{u}^{(z\mu)}\_{ij}(\bar{r}) = -\frac{ez\_i\left(\mu\_j\cdot \hat{r}\right) }{|\bar{r}|^2} \left( \mathcal{S}(q) - q\mathcal{S}^{\prime}(q) \right)
$$

where $\hat{r} = \bar{r}/|\bar{r}|$, and the dipole-dipole interaction by

$$
\tilde{u}^{\mu\mu}\_{ij}(\bar{r}) = -\left ( \frac{3 ( \boldsymbol{\mu}\_i \cdot \hat{r} ) \left(\boldsymbol{\mu}\_j\cdot\hat{r}\right) - \boldsymbol{\mu}\_i\cdot\boldsymbol{\mu}\_j }{|\bar{r}|^3}\right) \left( \mathcal{S}(q) - q\mathcal{S}^{\prime}(q)  + \frac{q^2}{3}\mathcal{S}^{\prime\prime}(q) \right) - \frac{\left(\boldsymbol{\mu}\_i\cdot\boldsymbol{\mu}\_j\right)}{|\bar{r}|^3}\frac{q^2}{3}\mathcal{S}^{\prime\prime}(q).
$$

**Warning:**

 - The `zahn` and `fennell` approaches have undefined dipolar self-energies (see next section) and are therefore not recommended for dipolar systems. 


### Self-energies

When using `coulomb` or `multipole`, an electrostatic self-energy term is automatically
added to the Hamiltonian. The monopole and dipole contributions are evaluated according to

$$
U\_{self} = -\frac{1}{2}\sum\_i^N\sum\_{\ast\in\{z,\mu\}} \lim\_{|\bar{r}\_{ii}|\to 0}\left( u^{(\ast\ast)}\_{ii}(\bar{r}\_{ii})
- \tilde{u}^{(\ast\ast)}\_{ii}(\bar{r}\_{ii}) \right )
$$

where no tilde indicates that $\mathcal{S}(q)\equiv 1$ for any $q$.


### Ewald Summation

If type is `ewald`, terms from reciprocal space and surface energies are automatically added (in addition to the previously mentioned self- and real space-energy) to the Hamiltonian which activates the additional keywords:

`type=ewald`          | Description
--------------------- | ---------------------------------------------------------------------
`ncutoff`             | Reciprocal-space cutoff (unitless)
`epss=0`              | Dielectric constant of surroundings, $\varepsilon_{surf}$ (0=tinfoil)
`ewaldscheme=PBC`     | Periodic (`PBC`) or isotropic periodic ([`IPBC`](http://doi.org/css8)) boundary conditions
`spherical_sum=true`  | Spherical/ellipsoidal summation in reciprocal space; cubic if `false`.
`debyelength=`$\infty$| Debye length (Å)

The added energy terms are:

$$
U_{\text{reciprocal}} = \frac{2\pi f}{V} \sum_{ {\bf k} \ne {\bf 0}} A\_k \vert Q^{q\mu} \vert^2 
$$

$$
U_{\text{surface}} = \frac{1}{4\pi\varepsilon_0\varepsilon_r}\frac{ 2\pi }{ (2\varepsilon_{surf} + 1) V }
\left(
\left|\sum_{j}q\_j\bar{r}\_j\right|^2 + 2 \sum\_j q\_i \bar{r}\_j \cdot \sum\_j \boldsymbol{\mu}\_j + \left| \sum\_j \boldsymbol{\mu}\_j \right|^2
\right )
$$

where

$$
    f = \frac{1}{4\pi\varepsilon_0\varepsilon_r} \quad\quad V=L_xL_yL_z
$$

$$
A\_k = \frac{e^{-( k^2 + \kappa^2 )/4\alpha^2}}{k^2}
\quad \quad Q^{q\mu} = Q^{q} + Q^{\mu}
$$

$$
Q^{q} = \sum_{j}q\_je^{i({\bf k}\cdot {\bf r}\_j)} \quad Q^{\mu} = \sum_{j}i({\boldsymbol{\mu}}\_j\cdot {\bf k})  e^{i({\bf k}\cdot {\bf r}\_j)}
$$

$$
\bar{k} = 2\pi\left( \frac{n_x}{L_x} , \frac{n_y}{L_y} ,\frac{n_z}{L_z} \right)\quad \bar{n} \in \mathbb{Z}^3
$$

Like many other electrostatic methods, the Ewald scheme also adds a self-energy term as described above.
In the case of isotropic periodic boundaries (`ipbc=true`), the orientational degeneracy of the
periodic unit cell is exploited to mimic an isotropic environment, reducing the number
of wave-vectors to one fourth compared with 3D PBC Ewald.
For point charges, [IPBC](http://doi.org/css8) introduce the modification,

$$
Q^q = \sum\_j q\_j \prod\_{\alpha\in\{x,y,z\}} \cos \left( \frac{2\pi}{L\_{\alpha}} n\_{\alpha} r\_{\alpha,j} \right)
$$

while for point dipoles (currently unavailable),

$$
Q^{\mu} = \sum\_j \bar{\mu}\_j
\cdot \nabla\_j
\left( \prod\_{ \alpha \in \{ x,y,z \} } \cos \left ( \frac{2\pi}{L\_{\alpha}} n\_{\alpha} \bar{r}\_{\alpha,j} \right ) \right )
$$

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
`dz=0.2`          | $z$ resolution (Å)
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

where $a\_j$ is the radius of the non-polar particle and $\alpha\_j$ is set in
the atom topology, `alphax`.
For non-polar particles in a polar medium, $\alpha\_i$ is a negative number.
For more information, see
[J. Israelachvili's book, Chapter 5.](https://www.sciencedirect.com/science/book/9780123751829)

`ionalpha`           | Description
-------------------- | ---------------------------------------
`epsr`               | Relative dielectric constant of medium

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
`rc`   | Width, $r\_c$ (Å)
`wc`   | Decay range, $w\_c$ (Å)

### Assorted Short Ranged Potentials

The potentials below are often used to keep particles apart and/or to introduce stickiness.
The atomic interaction parameters, e.g., $\sigma\_i$ and $\epsilon\_i$, are taken from the
topology.

Type             | Atomic parameters | $u(r)$ (non-zero part)
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

Rule                  | Description        | Formula
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
      - Na Cl: {eps: 0.2, sigma: 2}
      - K Cl: { ... }
- hertz:
    mixing: LB
    eps: eps_hz
    custom:
      - Na Cl: {eps_hz: 0.2, sigma: 2}
- hardsphere:
    mixing: arithmetic
    sigma: sigma_hs
    custom:
      - Na Cl: {sigma_hs: 2}
~~~

### SASA (pair potential)

This calculates the surface area of two intersecting particles or radii $R$ and $r$
to estimate an energy based on transfer-free-energies (TFE) and surface tension.
The total surface area is calculated as

$$
    A = 4\pi \left ( R^2 + r^2 \right ) - 2\pi \left (  Rh\_1 + rh\_2  \right )
$$

where $h\_1$ and $h\_2$ are the heights of the spherical caps comprising the lens
formed by the overlapping spheres. For complete overlap, or when far apart, the
full area of the bigger sphere or the sum of both spheres are returned.
The pair-energy is calculated as:

$$
    u_{ij} = A \left ( \gamma\_{ij} + c\_s \varepsilon\_{\text{tfe},ij} \right )
$$

where $\gamma\_{ij}$ and $\varepsilon\_{\text{tfe},ij}$ are the arithmetic means of
 `tension` and `tfe` provided in the atomlist.

Note that SASA is strictly not additive and this pair-potential is merely
a poor-mans way of approximately taking into account ion-specificity and
hydrophobic/hydrophilic interactions. Faunus offers also a full, albeit yet
experimental implementation of [Solvent Accessible Surface Area] energy.

`sasa`       | Description
------------ | ----------------------------------------------------------
`molarity`   | Molar concentration of co-solute, $c\_s$
`radius=1.4` | Probe radius for SASA calculation (Å)
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
`cutoff`     | Spherical cutoff distance

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
`e0`       | Vacuum permittivity [C²/J/m]
`inf`      | ∞ (infinity)
`kB`       | Boltzmann constant [J/K]
`kT`       | Boltzmann constant × temperature [J]
`Nav`      | Avogadro's number [1/mol]
`pi`       | π (pi)
`q1`,`q2`  | Particle charges [e]
`r`        | Particle-particle separation [Å]
`Rc`       | Spherical cut-off [Å]
`s1`,`s2`  | Particle sigma [Å]
`T`        | Temperature [K]

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
`e0`       | Vacuum permittivity [C²/J/m]
`inf`      | ∞ (infinity)
`kB`       | Boltzmann constant [J/K]
`kT`       | Boltzmann constant × temperature [J]
`Nav`      | Avogadro's number [1/mol]
`pi`       | π (pi)
`q`        | Particle charge [e]
`s`        | Particle sigma [Å]
`x`,`y`,`z`| Particle positions [Å]
`T`        | Temperature [K]

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

### Gouy Chapman

By setting `function=gouychapman`, an electric potential from a uniformly, charged plane
in a 1:1 salt solution is added; see _e.g._ the book _Colloidal Domain_ by Evans and Wennerström, 1999.
If a surface potential, $\varphi\_0$ is specified,

$$
\rho = \sqrt{\frac{2 c\_0}{\pi \lambda\_B} } \sinh ( \beta e \varphi\_0 / 2 )
$$
while if instead a surface charge density, $\rho$, is given,
$$
\beta e \varphi\_0 = 2\mbox{asinh} \left ( \rho \sqrt{\frac{\pi \lambda\_B} {2 c\_0}} \right )
$$
where $\lambda\_B$ is the Bjerrum length. With $\Gamma\_0 = \tanh{ \beta e \varphi\_0 / 4 }$
the final, non-linearized external potential is:
$$
\beta e \phi\_i = 2 \ln \left ( \frac{1+\Gamma\_0e^{-\kappa r\_{z,i}}}{1-\Gamma\_0 e^{-\kappa r\_{z,i}}} \right )
$$
where
$z\_i$ is the particle charge;
$e$ is the electron unit charge;
$\kappa$ is the inverse Debye length;
and $r\_{z,i}$ is the distance from the charged $xy$-plane which is always placed at the minimum
$z$-value of the simulation container (normally a slit geometry).
Fluctuations of the simulation cell dimensions are respected.

The following parameters should be given under `constants`;
the keywords `rho`, `rhoinv`, and `phi0` are mutually exclusive.

`constants`        | Description
-----------------  | ----------------------------------------------------------------------
`molarity`         | Molar 1:1 salt concentration (mol/l)
`epsr`             | Relative dielectric constant
`rho`              | Charge per area (1/eÅ²)
`rhoinv`           | Area per charge (eÅ²) if `rho` nor `phi0` are given
`phi0`             | Unitless surface potential, $\beta e \varphi\_0$, if `rho` or `rhoinv` not given
`linearise=false`  | Use linearised Poisson-Boltzmann approximation?


## Custom Group-Group Potential

For two different or equal molecule types (`name1`, `name2`), this adds a user-defined energy function
given at run-time. In addition to physical constants, the following varibles are available:

`constants`        | Description
-----------------  | --------------------------------
`R`                | Mass center separation (Å)
`Z1`               | Average net-charge of group 1
`Z2`               | Average net-charge of group 2

When used together with regular non-bonded interactions, this can e.g. be used to bias simulations.
In the following example, we _subtract_ a Yukawa potential and the bias can later be removed by
re-weighting using information from the `systemenergy` analysis output.

~~~ yaml
custom-groupgroup:
    name1: charged_colloid
    name2: charged_colloid
    constants: { bjerrum_length: 7.1, debye_length: 50 }
    function: >
        -bjerrum_length * Z1 * Z2 / R * exp(-R / debye_length)
~~~


## Bonded Interactions

Bonds and angular potentials are added via the keyword `bondlist` either directly
in a molecule definition (topology) for intra-molecular bonds, or in `energy->bonded`
where the latter can be used to add inter-molecular bonds:

~~~ yaml
moleculelist:
    - water: # TIP3P
        structure: "water.xyz"
        bondlist: # index relative to molecule
            - harmonic: { index: [0,1], k: 5024, req: 0.9572 }
            - harmonic: { index: [0,2], k: 5024, req: 0.9572 }
            - harmonic_torsion: { index: [1,0,2], k: 628, aeq: 104.52 }
energy:
    - bonded:
        bondlist: # absolute index; can be between molecules
           - harmonic: { index: [56,921], k: 10, req: 15 }
~~~

$\mu V T$ ensembles and Widom insertion are currently unsupported for molecules with bonds.

The following shows the possible bonded potential types:

### Harmonic

`harmonic`     | Harmonic bond
-------------- | -------------------------------------------
`k`            | Harmonic spring constant (kJ/mol/Å²)
`req`          | Equilibrium distance (Å)
`index`        | Array with _exactly two_ indices (relative to molecule)

$$
u(r) = \frac{1}{2}k(r-r\_{\mathrm{eq}})^2
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
-\frac{1}{2} k r_{\mathrm{max}}^2 \ln \left [ 1-(r/r_{\mathrm{max}})^2 \right ],  & \text{if } r < r_{\mathrm{max}} \\
\infty, & \text{if } r \geq r_{\mathrm{max}}
\end{cases}
$$

It is recommended to only use the potential if the initial configuration is near equilibrium, which prevalently depends on the value of `rmax`.
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
with the short ranged Weeks-Chandler-Andersen (wca) repulsive potential. This potential is particularly useful in combination with the `nonbonded_cached` energy.

$$
     u(r) =
  \begin{cases} 
   -\frac{1}{2} k r_{\mathrm{max}}^2 \ln \left [ 1-(r/r_{\mathrm{max}})^2 \right ] + u_{\mathrm{wca}}, & \text{if } 0 < r \leq 2^{1/6}\sigma \\
   -\frac{1}{2} k r_{\mathrm{max}}^2 \ln \left [ 1-(r/r_{\mathrm{max}})^2 \right ],       & \text{if } 2^{1/6}\sigma < r < r_{\mathrm{max}} \\
   \infty, & \text{if } r \geq r_{\mathrm{max}}
  \end{cases}
$$

It is recommended to only use this potential if the initial configuration is near equilibrium, which prevalently depends on the value of `rmax`.
Should one insist on conducting simulations far from equilibrium, a large displacement parameter is recommended to reach finite energies.

### Harmonic torsion

`harmonic_torsion` | Harmonic torsion
------------------ | -------------------------------------------
`k`                | Harmonic spring constant (kJ/mol/rad²)
`aeq`              | Equilibrium angle $\alpha_{\mathrm{eq}}$ (deg)
`index`            | Array with _exactly three_ indices (relative to molecule)

$$
u(\alpha) = \frac{1}{2}k(\alpha - \alpha_{\mathrm{eq}})^2
$$
where $\alpha$ is the angle between vector 1→0 and 1→2 (numbers refer to the position in `index`).

### Cosine based torsion (GROMOS-96)

`gromos_torsion`   | Cosine based torsion
------------------ | -------------------------------------------
`k`                | Force constant (kJ/mol)
`aeq`              | Equilibrium angle $\alpha_{{\mathrm{eq}}}$ (deg)
`index`            | Array with _exactly three_ indices (relative to molecule)

$$
u(\alpha) = \frac{1}{2}k(\cos(\alpha) - \cos(\alpha_{{\mathrm{eq}}}))^2
$$
where $\alpha$ is the angle between vector 1→0 and 1→2 (numbers refer to the position in `index`).


### Proper periodic dihedral

`periodic_dihedral` | Proper periodic dihedral
------------------- | -------------------------------------------
`k`                 | Force constant (kJ/mol)
`n`                 | Periodicity (multiplicity) of the dihedral (integer)
`phi`               | Phase angle $\phi_{\mathrm{syn}}$ (deg)
`index`             | Array with _exactly four_ indices (relative to molecule)

$$
u(r) = k(1 + \cos(n\phi - \phi_{\mathrm{syn}}))
$$
where $\phi$ is the angle between the planes constructed from the points 0,1,2 and 1,2,3 (numbers refer to the position in `index`).


### Improper harmonic dihedral

`harmonic_dihedral` | Improper harmonic dihedral
------------------- | -------------------------------------------
`k`                 | Harmonic spring constant (kJ/mol/rad²)
`deq`               | Equilibrium angle $\alpha_{\mathrm{eq}}$ (deg)
`index`             | Array with _exactly four_ indices (relative to molecule)

$$
u(\phi) = \frac{1}{2}k(\phi - \phi_{\mathrm{eq}})^2
$$
where $\phi$ is the angle between the planes constructed from the points 0,1,2 and 1,2,3 (numbers refer to the position in `index`).


## Geometrical Confinement

`confine`    | Confine molecules to a sub-region
------------ | -------------------------------------------
`type`       | Confinement geometry: `sphere`, `cylinder`, or `cuboid`
`molecules`  | List of molecules to confine (names)
`com=false`  | Apply to molecular mass center
`k`          | Harmonic spring constant in kJ/mol or `inf` for infinity

Confines `molecules` in a given region of the simulation container by applying a harmonic potential on
exterior atom positions, $\mathbf{r}\_i$:

$$
U = \frac{1}{2} k \sum\_{i}^{\text{exterior}} f\_i
$$

where $f\_i$ is a function that depends on the confinement `type`,
and $k$ is a spring constant. The latter
may be _infinite_ which renders the exterior region strictly inaccessible.
During equilibration it is advised to use a _finite_ spring constant to drive exterior particles inside the region.
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

`sasa`       | SASA Transfer Free Energy
------------ | --------------------------------------------
`radius=1.4` | Probe radius for SASA calculation (Å)
`molarity`   | Molar concentration of co-solute
`dense=true` | Flag specifying if a dense or a sparse version of a cell list container is used
`slices=25`  | Number of slices per particle when calculating SASA (the more, the more precise)

Calculates the free energy contribution due to

1. atomic surface tension
2. co-solute concentration (typically electrolytes)

via a [SASA calculation](http://dx.doi.org/10/dbjh) for each particle.
The energy term is:

$$
    U = \sum_i^N A_{\text{sasa},i} \left ( \gamma_i + c_s \varepsilon_{\text{tfe},i} \right )
$$

where $c\_s$ is the molar concentration of the co-solute;
$\gamma\_i$ is the atomic surface tension; and $\varepsilon_{\text{tfe},i}$ the atomic transfer free energy,
both specified in the atom topology with `tension` and `tfe`, respectively.
Will use cell lists if a geometry is either `cuboid` or `sphere`.
The `dense` option specifies if a dense implementation 
(memory heavy but faster) or a sparse one (slightly slower but light) of a cell list container will be used.

### Alternative schemes

Scheme           | PBC | Cell list | Note
-----------------|:---:|:---------:|-----------------------------
`sasa`           |  ✓  |     ✓     | Default
`sasa_reference` |  ✓  |     ✓     | Use for debugging
`freesasa`       |     |           | Uses [FreeSASA](https://freesasa.github.io/)


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
`f0`             |  Penalty energy increment (kT)
`update`         |  Interval between scaling of `f0`
`scale`          |  Scaling factor for `f0`
`nodrift=true`   |  Suppress energy drift
`quiet=false`    |  Set to true to get verbose output
`file`           |  Name of saved/loaded penalty function
`overwrite=true` |  If `false`, don't save final penalty function
`histogram`      |  Name of saved histogram (optional)
`coords`         |  Array of _one or two_ coordinates

The coordinate, $\mathcal{X}$, can be freely composed by one or two
of the types listed in the next section (via `coords`).


### Reaction Coordinates

The following reaction coordinates can be used for penalising the energy and can further
be used when analysing the system (see Analysis).
Please notice that atom id's are determined by the order of appearance in the `atomlist` whereas molecular id's
follow the order of insertion specified in `insertmolecules`.

General keywords  | Description 
----------------- | -------------------------------------------------------------------
`index`           | Atom index, atom id or group index
`indexes`         | Array of atomic indexes (`[a,b]` or `[a,b,c,d]`)
`range`           | Array w. [min:max] value
`resolution`      | Resolution along the coordinate (Å)
`dir`             | Axes of the reaction coordinate, $e.g.$, `[1,1,0]` for the $xy$-plane  

#### Atom Properties

`coords=[atom]`| Property
-------------- | --------------------------------------------------------------------------
`x`, `y` or `z`| $x$-, $y$- or $z$-coordinate of the $i$th particle, $i$=`index`
`q`            | Charge of the $i$th particle, $i$=`index`
`R`            | Distance of the $i$th particle from the center of the simulation cell, $i$=`index`
`N`            | Number of atoms of id=`index`

#### Molecule Properties

`coords=[molecule]`         | Property
--------------------------- | ---------------------------------------------------------------------------
`active`                    | If molecule is active (1) or inactive (0); for GCMC ensembles
`angle`                     | Angle between instantaneous principal axis and given `dir` vector
`com_x`, `com_y` or `com_z` | Mass-center coordinates
`confid`                    | Conformation id corresponding to frame in `traj` (see molecular topology).
`end2end`                   | Distance between first and last atom
`Rg`                        | Radius of gyration
`mu_x`, `mu_y` or `mu_z`    | Molecular dipole moment components
`mu`                        | Molecular dipole moment scalar ($e$Å/charge)
`muangle`                   | Angle between dipole moment and given `dir` vector
`N`                         | Number of atoms in group
`Q`                         | Monopole moment (net charge)
`atomatom`                  | Distance along `dir` between 2 atoms specified by the `indexes` array
`cmcm`                      | Absolute mass-center separation between group indexes `a` and `b` or atomic indexes `a`–`b` and `c`–`d`
`cmcm_z`                    | $z$-component of `cmcm`
`mindist`                   | Minimum distance between particles of id `indexes[0]` and `indexes[1]`
`L/R`                       | Ratio between height and radius of a cylindrical vesicle
`Rinner`                    | Average $d$ of id=`indexes[0]` for particles having a smaller $d$ than id=`indexes[1]`

Notes:

- the molecular dipole moment is defined with respect to the mass-center
- for `angle`, the principal axis is the eigenvector corresponding to the smallest eigenvalue of the gyration tensor
- `Rinner` can be used to calculate the inner radius of cylindrical or spherical vesicles. $d^2=\bf{r} \cdot$`dir` where
$\bf{r}$ is the position vector
- `L/R` can be used to calculate the bending modulus of a cylindrical lipid vesicle
- `Rg` is calculated as the square-root of the sum of the eigenvalues of the gyration tensor, $S$. 
$$
S = \frac{1}{\sum_{i=1}^{N} m_{i}} \sum_{i=1}^{N} m_{i} \bf{t_i} \bf{t_i^T}
$$
where $\bf{t_i} = \bf{r_i} - \bf{cm}$, $\bf{r_i}$ is the coordinate of the $i$th atom, $m_i$ is the mass of the $i$th atom, $\bf{cm}$ is the
mass center of the group and $N$ is the number of atoms in the molecule.

#### System Properties

`coords=[system]`         | Property
------------------------- | -----------------------------------------------
`V`                       | System volume
`Q`                       | System net-charge
`Lx`, `Ly` or `Lz`        | Side length of the cuboidal simulation cell
`height`                  | Alias for `Lz`
`radius`                  | Radius of spherical or cylindrical geometries
`N`                       | Number of active particles
`mu`                      | System dipole moment scalar (𝑒Å)
`mu_x`, `mu_y` or `mu_z`  | System dipole moment components (𝑒Å)

The enclosing cuboid is the smallest cuboid that can contain the geometry.
For example, for a cylindrical simulation container, `Lz` is the height
and `Lx=Ly` is the diameter.


### Multiple Walkers with MPI

If compiled with MPI, the master process collects the bias function from all nodes
upon penalty function `update`.
The _average_ is then re-distributed, offering [linear parallelization](http://dx.doi.org/10/b5pc4m)
of the free energy sampling. It is crucial that the walk in coordinate space differs in the different
processes, e.g., by specifying a different random number seed; start configuration; or displacement parameter.
File output and input are prefixed with `mpi{rank}`.

The following starts all MPI processes with the same input file, and the MPI prefix is automatically
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

