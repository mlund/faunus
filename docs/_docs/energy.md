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
    - isobaric: { P/atm: 1 }
    - sasa: { molarity: 0.2, radius: 1.4 }
    - confine: { type: sphere, radius: 10,
                 molecules: ["water"] }
    - nonbonded_coulomblj:
        lennardjones: { mixing: LB }
        coulomb: { type: plain, epsr: 80, cutoff: 10 }
    - ...
~~~

**Note:**
_Energies_ in MC may contain implicit degrees of freedom, _i.e._ be temperature-dependent,
effective potentials. This is inconsequential for sampling
density of states, but care should be taken when interpreting derived functions such as
energies, entropies, pressure etc.
{: .notice--info}

## External Pressure <a name="isobaric"></a>

`isobaric`   | Description
-------------|-------------------------------------------------------
`P/unit`     | External pressure where unit can be `mM`, `atm`, or `Pa`.

This adds the following pressure term[^frenkel] to the Hamiltonian, appropriate for
MC moves in $\ln V$:

$$
U = PV - k_BT\left ( N + 1 \right ) \ln V
$$

where $N$ is the total number of molecules and atomic species.

[^frenkel]: _Frenkel and Smith, 2nd Ed., Chapter 5.4_.

## Nonbonded Interactions <a name="nonbonded"></a>

`energy`               | $u_{ij}$
---------------------- | --------------------------
`nonbonded_coulomblj`  | `coulomb`+`lennardjones`
`nonbonded_coulombhs`  | `coulomb`+`hardsphere`
`nonbonded_coulombwca` | `coulomb`+`wca`

This term loops over pairs of atoms, $i$, and $j$, summing a given pair-wise additive potential, $u_{ij}$,

$$ U = \sum_{i=0}^{N-1}\sum_{j=i+1}^N u_{ij}(\textbf{r}_j-\textbf{r}_i)$$

**Note:** the pair-potential can be a combination of several potentials, but this must currently be defined at _compile
time_ and cannot be arbitrarily selected from input. It is straight forward to add more by editing the class
`Energy::Hamiltonian` found in `src/energy.h` and then re-compile.
{: .notice--info}

Below is a description of possible pair-potentials and their configuration.

### Coulomb

 `coulomb`   |  Description
 ----------- |  -------------------------------------------------
 `type`      |  Coulomb type, see below
 `cutoff`    |  Spherical cutoff, $R_c$ after which the potential is zero
 `epsr`      |  Relative dielectric constant of the medium

This is a multipurpose pair potential that handles several variants of electrostatic
pair-potentials. Beyond a spherical cutoff, $R_c$, the potential is zero while if
below,

$$
u_{ij} = \frac{e^2 z_i z_j }{ 4\pi\epsilon_0\epsilon_r r_{ij} }\mathcal{S}(q)
$$

where $\mathcal{S}(q=r/R_c)$ is a splitting function:

 `type`         | $\mathcal{S}(q)$                                    | Keywords      | Ref.
 -------------- | --------------------------------------------------- | ------------- | ----------------------
 `none`         | $0$                                                 |               |
 `plain`        | $1$                                                 |               | [doi](http://doi.org/ctnnsj)
 `ewald`        | $\text{erfc}(\alpha R_cq)$                          | `alpha`       | [doi](http://dx.doi.org/10.1063/1.481216)
 `wolf`         | $\text{erfc}(\alpha R_cq)-\text{erfc}(\alpha R_c)q$ | `alpha`       | [doi](http://doi.org/cfcxdk)
 `yukawa`       | $e^{-\kappa R_c q}-e^{-\kappa R_c}$                 | `debyelength` | [ISBN](https://isbnsearch.org/isbn/0486652424)
 `yonezawa`     | $1+\text{erfc}(\alpha R_c)q+q^2$                    | `alpha`       | [doi](http://dx.doi.org/10/j97)
 `qpotential`   | $\prod_{n=1}^{order}(1-q^n)$                        | `order=300`   | [ISBN](http://goo.gl/hynRTS) (Paper V)
 `fanourgakis`  | $1-\frac{7}{4}q+\frac{21}{4}q^5-7q^6+\frac{5}{2}q^7$|               | [doi](http://doi.org/f639q5)
 `fennel`       | $\scriptstyle \text{erfc}(\alpha R_c q)-\text{erfc}(\alpha R_c)q + ( q -1 ) q \left( \text{erfc}(\alpha R_c) + \frac{2\alpha R_c}{\sqrt{\pi}} e^{-\alpha^2 R_c^2} \right)$ | `alpha`| [doi](http://doi.org/bqgmv2)
 `reactionfield`| $1+\frac{\epsilon_{RF}-\epsilon_{r}}{2\epsilon_{RF}+\epsilon_{r}}q^3-3\frac{\epsilon_{RF}}{2\epsilon_{RF}+\epsilon_{r}}q$ | `epsrf` | [doi](http://doi.org/dbs99w)

**Note:** Internally $\mathcal{S}(q)$ is _splined_ whereby all types evaluate at similar speed.
{: .notice--info}

#### Ewald Summation

If type is `ewald`, terms from reciprocal space; surface energies; and
self energies are automatically added to the Hamiltonian, activating additional keywords:

`type=ewald`         | Description
-------------------- | --------------------------------------
`cutoffK`            | Reciprocal-space cutoff
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
\quad\quad
Q^{q\mu} = \sum_{j}q_j + i({\boldsymbol{\mu}}_j\cdot {\bf k})  e^{i({\bf k}\cdot {\bf r}_j)}
$$

$$
{\bf k} = 2\pi\left( \frac{n_x}{L_x} , \frac{n_y}{L_y} ,\frac{n_z}{L_z} \right),\;\; {\bf n} \in \mathbb{Z}^3
$$


### Hard Sphere
`hardsphere`

The hard sphere potential does not take any input. Radii are read from the atomlist at the beginning of the simulation.


### Lennard-Jones

`lennardjones` |  Description
-------------  |  ------------------------------------------------
`mixing=LB`    |  Mixing rule. `LB`
`ljcustom`     |  Custom $\epsilon$ and $\sigma$ combinations

The Lennard-Jones potential has the form:

$$
u_{ij}^{\text{LJ}} = 4\epsilon_{ij} \left (
    \left ( \frac{\sigma_{ij}} {r_{ij})} \right )^{12} - \left ( \frac{\sigma_{ij}}{r_{ij})}\right )^6 \right )
$$

where the default mixing rule is Lorentz-Berthelot (`LB`):

$$
\sigma_{ij} = \frac{\sigma_i+\sigma_j}{2} \quad \textrm{and} \quad
\epsilon_{ij} = \sqrt{\epsilon_i \epsilon_j}
$$

### Weeks-Chandler-Andersen

`wca`          |  Description
-------------  |  ------------------------------------------------
`mixing=LB`    |  Mixing rule; only `LB` available.
`ljcustom`     |  Custom $\epsilon$ and $\sigma$ combinations

Like Lennard-Jones but cut and shifted to zero at the minimum, forming a
purely repulsive potential [[doi](http://dx.doi.org/ct4kh9)],

$$
u_{ij} = u_{ij}^{\text{LJ}} + \epsilon_{ij} \quad \textrm{for} \quad r<2^{1/6}\sigma_{ij}
$$

## Bonded Interactions

Bonds and angular potentials are added via the keyword `bondlist` either directly
in a molecule definition (topology) or in energy/bonded where the latter can be
used to add inter-molecular bonds:

~~~ yaml
moleculelist:
    - water:
        structure: water.xyz
        bondlist: # index relative to molecule
            - harmonic: { index: [0,1], k: 100, req: 1.0 }  
            - harmonic: { index: [0,2], k: 100, req: 1.0 }
energy:
    - bonded:
        bondlist: # absolute index
           - harmonic: { index: [56,921], k: 10, req: 15 }
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
`index`        | array with _exactly two_ index (relative to molecule)

$$
u(r) = \frac{1}{2}k(r-r_{eq})^2
$$

### Finite Extensible Nonlinear Elastic

`fene`         | [Finite Extensible Nonlinear Elastic Potential](http://dx.doi.org/10.1103/PhysRevE.59.4248)
-------------- | -------------------------------------------
`k`            | Bond stiffness (kJ/mol/Å$^2$)
`rmax`         | Maximum separation, $r_m$ (Å)
`index`        | array with _exactly two_ index (relative to molecule)

$$
u(r) = -\frac{1}{2} k r_m^2 \ln \left [ 1-(r/r_m)^2 \right ]
$$

for $r < r_m$; infinity otherwise.

## Geometrical Confinement

`confine`    | Confine molecules to a sub-region
------------ | -------------------------------------------
`type`       | Confinement geometry: `sphere`, `cylinder`, or `cuboid`
`molecules`  | List of molecules to confine (names)
`k`          | Harmonic spring constant in kJ/mol or `inf` for infinity

Confines `molecules` in a given region of the simulation container by applying a harmonic potential on
exterior atom positions, $\mathbf{r}_i$:

$$
U = \frac{1}{2} k \sum^{exterior} f_i
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

via a [SASA calculation](http://dx.doi.org/10.1002/jcc.21844) for each atom.
The energy term is:

$$
U = \sum_i^N A_{\text{sasa},i} \left ( \gamma_i + c \cdot \varepsilon_{\text{tfe},i} \right )
$$

where $c$ is the molar concentration of the co-solute;
$\gamma_i$ is the atomic surface tension; and $\varepsilon_{\text{tfe},i}$ the atomic transfer free energy,
both specified in the atom topology with `tension` and `tfe`, respectively.


