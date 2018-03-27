---
---
<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
});
</script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>
[![Edit](https://img.shields.io/badge/Github-Improve_this_page-orange.svg)]({{site.github.repository_url}}/blob/master/docs/{{page.path}})

# Analysis

Faunus can perform on-the-fly analysis during simulation by allowing for an arbitrary number
of analysis functions to be added. The list of analysis is defined in the `analysis` section
at the top level input:

~~~ yaml
analysis:
    - systemenergy: {file: energy.dat, nstep: 500}
    - xtcfile: {file: traj.xtc, nstep: 1000}
    - widom:  {molecule: water, ninsert: 20, nstep: 50}
    - molrdf: {name1: water, name2: water, nstep: 100,
               dr: 0.1, dim: 3, file: rdf.dat}
    - ...
~~~

**Note:** all analysis methods support the `nstep` keyword that defines the interval between
sampling points. In addition all analysis provide output statistics of number of sample
points, and the relative run-time spent on the analysis.
{: .notice--info}

## Density

`density`   |  Description
----------- |  -------------------------------------------
`nstep=0`   |  Interval between samples

This calculates the average density, $\langle N_i/V \rangle$ of molecules and atoms
which may fluctuate in _e.g._ the isobaric ensemble or the Grand Canonical ensemble.
For atomic groups, densities of individual atom types are reported.

## Density Profile

`atomprofile`  | Description
-------------- | ---------------------------------------------
`nstep=0`      | Interval between samples
`atoms=[]`     | List of atom names to sample
`file`         | Output filename with profile
`dr=0.1`       | Radial resolution
`origo=[0,0,0]`| Center of the profile ($r=0$)

Calculates the summed density of `atoms` in spherical shells around
`origo` which by default is the center of the simulation box:

$$
\rho(r) = \frac{\langle N(r) \rangle}{4\pi r^2dr}
$$

This can be used to obtain charge profiles, measure excess pressure
etc.

## Structure

### Atomic $g(r)$

Samples the pair correlation function between atom id's _i_ and _j_,

$$
    g_{ij}(r) = \frac{ N_{ij}(r) }{ \sum_{r=0}^{\infty} N_{ij}(r) } \cdot \frac{ \langle V \rangle }{ V(r) }
$$

where $N_{ij}(r)$ is the number of observed pairs, accumulated over the
entire ensemble, in the separation
interval $[r, r+dr]$ and $V(r)$ is the corresponding volume element
which depends on dimensionality, `dim`.

`atomrdf`      |  Description
-------------- | ---------------------------------------------------------
`file`         |  Output file, two column
`name1`        |  Atom name 1
`name2`        |  Atom name 2
`dr=0.1`       |  $g(r)$ resolution
`dim=3`        |  Dimensions for volume element
`nstep=0`      |  Interval between samples

`dim` |  $V(r)$        
----- | ---------------
3     |  $4\pi r^2 dr$ 
2     |  $2\pi r dr$   
1     |  $dr$          

### Molecular $g(r)$

Same as `atomrdf` but for molecular mass-centers.

`molrdf`       |  Description
-------------- | ------------------------------
`file`         |  Output file, two column
`name1`        |  Molecule name 1
`name2`        |  Molecule name 2
`dr=0.1`       |  $g(r)$ resolution
`dim=3`        |  Dimensions for volume element
`nstep=0`      |  Interval between samples.

### Structure Factor

The isotropically averaged static structure factor between
$N$ point scatterers is calculated using the [Debye formula](http://doi.org/dmb9wm),

$$
    S(q) = 1 + \frac{2}{N} \left \langle
           \sum_{i=1}^{N-1}\sum_{j=i+1}^N \frac{\sin(qr_{ij})}{qr_{ij}}
           \right \rangle
$$

The selected `molecules` can be treated either as single
point scatterers (`com=true`) or as a group of individual
point scatterers of equal intensity, i.e. with a 
form factor of unity.

`scatter`   | Description
----------- | ------------------------------------------
`nstep`     | Interval with which to sample
`file`      | Output filename for $I(q)$
`molecules` | List of molecule names to sample (array)
`qmin`      | Minimum _q_ value (1/Å)
`qmax`      | Maximum _q_ value (1/Å)
`dq`        | _q_ spacing (1/Å)
`com=true`  | Treat molecular mass centers as single point scatterers


## Radius of Gyration

This calculate the radius of gyration, end-to-end distance, and related
fluctuations for all groups defined in `molecules`.

`polymershape`   | Description
---------------- | ----------------------------------------
`nstep`          | Interval with which to sample
`molecules`      | List of molecule names to sample (array)


## Molecular Multipoles

Calculates average molecular multipolar moments and their fluctuations.

`multipole`    | Description
-------------- | ----------------------
`nstep`        | Interval between samples.

## Electric Multipole Distribution

This will analyse the electrostatic energy between two groups as
a function of their mass center separation. Sampling consists of
the following:

1. The exact electrostatic energy is calculated by explicitly summing
   Coulomb interactions between charged particles
2. Each group - assumed to be a molecule - is translated into a
   multipole (monopole, dipole, quadrupole)
3. Multipolar interaction energies are calculated, summed, and tabulated
   together with the exact electrostatic interaction energy. Ideally
   (infinite number of terms) the multipoles should capture full
   electrostatics

The points 1-3 above will be done as a function of group-to-group
mass center separation, $R$ and moments
on molecule $a$ and $b$ with charges $q_i$ in position $\boldsymbol{r}_i$
with respect to the mass center are calculated according to:



$$
    q_{a/b} = \sum_i q_i \quad \quad \boldsymbol{\mu}_{a/b} = \sum_i q_i\mathbf{r_i}
$$

$$
    \boldsymbol{Q}_{a/b} = \frac{1}{2} \sum_i q_i\mathbf{r_i} \mathbf{r_i}^T
$$

And, omitting prefactors here, the energy between molecule $a$ and $b$ at $R$ is:

$$
    u_{\text{ion-ion}} = \frac{q_aq_b}{R} \quad \quad
    u_{\text{ion-dip}} = \frac{q_a \boldsymbol{\mu}_b \boldsymbol{R}}{R^3} + ...
$$

$$
\begin{aligned}
    u_{\text{dip-dip}}  = & \frac{\boldsymbol{\mu_a}\boldsymbol{\mu_b}   }{ R^3  }
        - \frac{3 (\boldsymbol{\mu_a} \cdot \boldsymbol{R}) ( \boldsymbol{\mu_b}\cdot\boldsymbol{R})  }{R^5}   \\
    u_{\text{ion-quad}} = & \frac{ q_a \boldsymbol{R}^T \boldsymbol{Q}_b \boldsymbol{R} }{R^5}
        - \frac{  q_a \mbox{tr}(\boldsymbol{Q}_b) }{R^3} + ... \\
    u_{\text{total}}    = & u_{\text{ion-ion}} + u_{\text{ion-dip}} + u_{\text{dip-dip}} + u_{\text{ion-quad}}\\
    u_{\text{exact}}    = & \sum_i^a\sum_j^b \frac{q_iq_j}{ | \boldsymbol{r_i} - \boldsymbol{r_j}  |    }
\end{aligned}
$$

During simulation, the above terms are thermally averaged over angles, co-solute degrees of freedom etc.
Note also that the moments are defined with respect to the *mass* center, not *charge* center.
While for most macromolecules there is only a minor difference between the two,
the latter is more appropriate and is planned for a future update.
A simply way to mimic this, is to assign zero mass to all neutral
atoms in the molecule.

The input keywords are:

`multipoledist`  | Description
---------------- | ----------------------------------------------------
`nstep`          | Interval between samples
`file`           | Output file name
`molecules`      | Array with exactly two molecule names, $a$ and $b$
`dr=0.2`         | Distance resolution (Å) along $R$.
 

## Save State

`savestate`    |  Description
-------------- | ---------------------------------------------------------
`file`         |  File to save; format detected by file extension: `pqr`, `aam`, `gro`, `xyz`, `state`
`nstep=-1`     |  Interval between samples. If -1, save at end of simulation

Saves the current configuration or the system state to file.
If the suffix is `state`, the system state is saved to a single
JSON file that can be used to restore the state.


## System Energy

Calculates the energy contributions from all terms in the Hamiltonian and
outputs to a file as a function of steps.
If filename ends with `.csv`, a comma separated value file will be saved,
otherwise a simple space separeted file with a single hash commented header line.
All units in $k_BT$.

`systemenergy`   |  Description
---------------- |  -------------------------------------------
`file`           |  Output filename for energy vs. step output
`nstep=0`        |  Interval between samples



## Virtual Volume Move

Performs a [virtual volume move](http://doi.org/cppxt6) by
scaling the simulation volume to $V+\Delta V$ along with
molecular mass centers and atomic positions. The excess pressure is evaluated
as a Widom average:

$$
    p^{ex} = \frac{k_BT}{\Delta V} \ln \left\langle e^{-\delta u / k_BT} \right\rangle_{NVT}
$$

For more advanced applications of volume perturbations - pressure tensors,
surface tension etc., see [here](http://doi.org/ckfh).

`virtualvolume` | Description
--------------- | -------------------------------------
`dV`            | Volume perturbation (angstrom cubed)
`nstep`         | Interval between samples


## Widom Insertion

This will insert a non-perturbing ghost molecule into
the system and calculate a [Widom average](http://doi.org/dkv4s6)
to measure the free energy of the insertion process, _i.e._ the
excess chemical potential:

$$
    \mu^{ex} = -k_BT \ln \left\langle e^{-\delta u/k_BT} \right\rangle_0
$$

where $\delta u$ is the energy change of the perturbation and the
average runs over the _unperturbed_ ensemble.
If the molecule has `atomic=true`, $\delta u$ includes the internal energy of the
inserted group. This is useful for example to calculate the mean excess activity
coefficient of a neutral salt pair.
Upon insertion, random positions and orientations are
generated. For use with rod-like particles on surfaces, the `absz`
keyword may be used to ensure orientations on only one
half-sphere.

**Important:**
One _inactive_ `molecule` must be added to the simulation using the `inactive`
keyword when inserting the initial molecules in the topology.
{: .notice--info}

`widom`       | Description
------------- | -----------------------------------------
`molecule`    | Name of _inactive_ molecule to insert (atomic or molecular)
`ninsert`     | Number of insertions per sample event
`dir=[1,1,1]` | Inserting directions
`absz=false`  | Apply `std::fabs` on all z-coordinates of inserted molecule
`nstep=0`      |  Interval between samples


## XTC trajectory

Generates a Gromacs XTC trajectory file with particle positions and box
dimensions as a function of steps.

`xtcfile`      |  Description
-------------- | ---------------------------------------------------------
`file`         |  Filename of output xtc file
`nstep=0`      |  Interval between samples.


