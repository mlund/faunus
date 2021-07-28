<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
});
</script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

# Analysis

Faunus can perform on-the-fly analysis during simulation by allowing for an arbitrary number
of analysis functions to be added. The list of analysis is defined in the `analysis` section
at the top level input:

~~~ yaml
analysis:
    - systemenergy: {file: energy.dat.gz, nstep: 500, nskip: 2000}
    - xtcfile: {file: traj.xtc, nstep: 1000}
    - widom:  {molecule: water, ninsert: 20, nstep: 50}
    - molrdf: {name1: water, name2: water, nstep: 100,
               dr: 0.1, dim: 3, file: rdf.dat}
    - ...
~~~

All analysis methods support the `nstep` keyword that defines the interval between
sampling points and the `nskip` keyword that defines the number of initial steps that are excluded from the analysis.
In addition all analysis provide output statistics of number of sample
points, and the relative run-time spent on the analysis.

## Density

### Bulk Density

`density`   |  Description
----------- |  -------------------------------------------
`nstep=0`   |  Interval between samples

This calculates the average density, $\langle N\_i/V \rangle$ of molecules and atoms
which may fluctuate in _e.g._ the isobaric ensemble or the Grand Canonical ensemble.
For atomic groups, densities of individual atom types are reported.
The analysis also files probability density distributions of atomic and polyatomic molecules
as well as of atoms involved in id transformations, _e.g._, acid-base equilibria.
The filename format is `rho-@name.dat`.

### Density Profile

`atomprofile`  | Description
-------------- | -----------------------------------------------------
`nstep=0`      | Interval between samples
`atoms=[]`     | List of atom names to sample; `[*]` selects all
`charge=false` | Calc. charge density instead of density
`file`         | Output filename with profile
`dr=0.1`       | Radial resolution
`origo=[0,0,0]`| Center of the profile ($r=0$)
`atomcom`      | Atom name; use the mass center of these atoms as origin
`dir=[1,1,1]`  | Direction along which the profile is calculated

Calculates the summed density of `atoms` in spherical, cylindrical or planar shells around
`origo` which by default is the center of the simulation box:

$$
\rho(r) = \frac{\langle N(r) \rangle}{V(r)}
$$

The sum of coefficients in `dir` determines the volume element normalisation:

$d_x+d_y+d_z$  | $V(r)$
-------------- | ----------------
3              | $4\pi r^2 dr$
2              | $2\pi r dr$
1              | $dr$

This can be used to obtain charge profiles, measure excess pressure
etc.

### Density Slice

`sliceddensity` | Description
--------------- | --------------------------------------------------------------
`atoms=[]`      | List of atom names to sample; `[*]` selects all 
`file`          | Output filename with profile
`dz=0.1`        | Resolution along _z_-axis
`atomcom`       | Atom name; use the mass center _z_ of these atoms as origin
`nstep`         | Interval between samples

Calculates the density in cuboidal slices of thickness _dz_ along the _z_ axis.
If an atom name is specified for the option `atomcom`, the _z_-position of each atom is calculated with respect to the center of mass of the atoms of the given type.

## Structure

### Atomic $g(r)$

Samples the pair correlation function between atom id's _i_ and _j_,

$$
    g_{ij}(r) = \frac{ N_{ij}(r) }{ \sum_{r=0}^{\infty} N_{ij}(r) } \cdot \frac{ \langle V \rangle }{ V(r) }
$$

where $N\_{ij}(r)$ is the number of observed pairs, accumulated over the
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
`slicedir`     |  Direction of the slice for quasi-2D/1D RDFs
`thickness`    |  Thickness of the slice for quasi-2D/1D RDFs

`dim` |  $V(r)$        
----- | ---------------
3     |  $4\pi r^2 dr$ 
2     |  $2\pi r dr$   
1     |  $dr$          

By specifying `slicedir`, the RDF is calculated only for atoms within a cylinder or slice of given `thickness`. For example, with `slicedir=[0,0,1]` and `thickness=1`, the RDF is calculated along _z_ for atoms within a cylinder of radius 1 Å. This quasi-1D RDF should be normalized with `dim=1`. Likewise, with `slicedir=[1,1,0]` and `thickness=2`, the RDF is calculated in the _xy_ plane for atoms with _z_ coordinates differing by less than 2 Å. This quasi-2D RDF should be normalized with `dim=2`.

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

### Dipole-dipole Correlation

Sample the dipole-dipole angular correlation, $\langle \pmb{\hat{\mu}}(0)\cdot \pmb{\hat{\mu}}(r) \rangle$, between dipolar atoms
and as a function of separation, _r_. In addition, the radial distribution function, $g(r)$ is sampled and saved to `{file}.gofr.dat`.

`atomdipdipcorr` |  Description
---------------- | ------------------------------
`file`           |  Output filename
`name1`          |  Atom name 1
`name2`          |  Atom name 2
`dr=0.1`         |  Angular correlation resolution
`dim=3`          |  Dimensions for volume element (affects only $g(r)$)
`nstep=0`        |  Interval between samples.


### Structure Factor

The isotropically averaged static structure factor between $N$ point scatterers is calculated using
the [Debye formula](http://doi.org/dmb9wm),

$$
    S(q) = 1 + \frac{2}{N} \left \langle
           \sum_{i=1}^{N-1}\sum_{j=i+1}^N \frac{\sin(qr_{ij})}{qr_{ij}}
           \right \rangle
$$

The selected `molecules` can be treated either as single point scatterers (`com=true`) or as a group of individual
point scatterers of equal intensity, i.e., with a  form factor of unity.

The computation of the structure factor is rather computationally intensive task, scaling quadratically with the number
of particles and linearly with the number of scattering vector mesh points. If OpenMP is available, multiple threads
may be utilized in parallel to speed it up the analysis.

`scatter`   | Description
----------- | ------------------------------------------
`nstep`     | Interval with which to sample
`file`      | Output filename for $I(q)$
`molecules` | List of molecule names to sample (array); `[*]` selects all 
`qmin`      | Minimum _q_ value (1/Å)
`qmax`      | Maximum _q_ value (1/Å)
`dq`        | _q_ spacing (1/Å)
`com=true`  | Treat molecular mass centers as single point scatterers
`pmax=15`   | Multiples of $(h,k,l)$ when using the `explicit` scheme
`scheme=explicit` | The following schemes are available: `debye`, `explicit`
`stepsave=false`  | Save every sample to disk

The `explicit` scheme is recommended for cuboids with PBC and the calculation is performed by explicitly averaging
the following equation over the 3+6+4 directions obtained by permuting the crystallographic index
`[100]`, `[110]`, `[111]` to define the scattering vector
$\mathbf{q} = 2\pi p/L(h,k,l)$ where $p=1,2,\dots,p\_{max}$.

$$
S(q) = \frac{1}{N} \left <
     \left ( \sum_i^N \sin(\mathbf{qr}\_i) \right )^2 +
     \left ( \sum_j^N \cos(\mathbf{qr}\_j) \right )^2
    \right >
$$

The sampled $q$-interval is always $\left [ 2\pi/L,\, 2\pi p\_{max} \sqrt{3} / L \right ]$,
$L$ being the box side length. Only cubic boxes have been tested, but the implementation respects cuboidal systems (untested).
For more information, see [doi:10.1063/1.449987](http://dx.doi.org/10.1063/1.449987).


### Atomic Inertia Eigenvalues

This calculates the inertia eigenvalues for all particles having a given id.
The inertia tensor is defined as

$$
I = \sum\_{i=1}^N m\_i ( \| \bf{t\_i} \|^2 \mathrm{I} - \bf{t\_i} \bf{t\_i}^T ) 
$$

where $\bf{t_i} = \bf{r_i} - \bf{cm}$, $\bf{r_i}$ is the coordinate of the $i$th particle, $\bf{cm}$ is the
position of the mass center of the whole group of atoms, $m_i$ is the molecular weight of the $i$th particle, 
$\bf{I}$ is the identity matrix and $N$ is the number of atoms.

`atominertia`       | Description
---------------- | ----------------------------------------
`nstep`          | Interval with which to sample
`index`          | Particle id

### Inertia Tensor

This calculates the inertia eigenvalues and the principal axis for a range of atoms within a molecular group of given index.
Atom coordinates are considered with respect to the mass center of the group. 
For protein complex, the analysis can be used to calculate the principal axes of the constituent monomers,
all originating at the mass center of the complex. 
The inertia tensor is defined as

$$
I = \sum_{i=1}^N m_i ( \| \bf{t_i} \|^2 \mathrm{I} - \bf{t_i} \bf{t_i}^T ) 
$$

where $\bf{t_i} = \bf{r_i} - \bf{cm}$, $\bf{r_i}$ is the coordinate of the $i$th particle, $\bf{cm}$ is the
position of the mass center of the whole group of atoms, $m_i$ is the molecular weight of the $i$th particle,
$\bf{I}$ is the identity matrix and $N$ is the number of atoms.

`inertia`    | Description
---------------- | ----------------------------------------
`nstep`          | Interval with which to sample
`indexes`        | Array defining a range of indexes within the molecule 
`index`          | Index of the molecular group

### Polymer Shape

This calculates the radius of gyration; end-to-end distance; and related
fluctuations for a molecular group. A histogram of the radius of gyration will
be saved to disk with the name `gyration_{molecule}.dat`.
The output nomenclature follows [IUPAC's recommendations](https://dx.doi.org/10/d6ff).
For further reading regarding gyration tensor analysis and shape, see:

- [doi:10.1021/ma00148a028](https://doi.org/10.1021/ma00148a028)
- [doi:10.1063/1.1730022](https://doi.org/10.1063/1.1730022)
- [wikipedia](https://en.wikipedia.org/wiki/Gyration_tensor)

From the principal moments, $\lambda$, the following shape descriptors are calculated:

- asphericity, $b = \lambda \_{{z}}^{{2}}-{\frac{1}{2}}\left(\lambda\_{{x}}^{{2}}+\lambda \_{{y}}^{{2}}\right)$.
- acylindricity, $c = \lambda \_{{y}}^{{2}}-\lambda \_{{x}}^{{2}}$
- relative shape anisotropy, $\kappa^2 ={\frac{3}{2}}{\frac  {\lambda\_{{x}}^{{4}}+\lambda\_{{y}}^{{4}}+\lambda\_{{z}}^{{4}}}{(\lambda \_{{x}}^{{2}}+\lambda \_{{y}}^{{2}}+\lambda \_{{z}}^{{2}})^{{2}}}}-{\frac{1}{2}}$



`polymershape`             | Description
-------------------------- | ----------------------------------------
`nstep`                    | Interval with which to sample
`molecule`                 | Molecule to sample
`histogram_resolution=0.2` | Rg resolution of histogram (Å)
`file`                     | Optionally save gyration tensor for each sample (.dat|.dat.gz)

Note: The ability to select several molecules (`molecules` keyword) was
removed in version 2.5. Instead, add multiple instances of `polymershape`.


### Molecular Conformation

For molecules that can have multiple conformations (using conformational swap moves), this
creates a histogram of observed conformations for a given molecule type.

`moleculeconformation` |  Description
---------------------- | -----------------
`nstep`                | Interval with which to sample
`molecule`             | Molecule name to sample


## Charge Properties

### Molecular Multipoles

Calculates average molecular multipolar moments and their fluctuations.

`multipole`    | Description
-------------- | ----------------------
`nstep`        | Interval between samples.

### Multipole Moments

For a range of atoms within a molecular group of given index, 
this calculates the total charge and dipole moment,
as well as the eigenvalues and the major axis of the quadrupole tensor.
Atom coordinates are considered with respect to the mass center of the group. 
For a protein complex, the analysis can be used to calculate, e.g., the dipole vectors of the constituent monomers,
all originating at the mass center of the complex. 
The quadrupole tensor is defined as

$$
Q = \frac{1}{2} \sum_{i=1}^N q_i ( 3 \bf{t_i} \bf{t_i}^T - \| \bf{t_i} \|^2 \mathrm{I}) 
$$

where $\bf{t_i} = \bf{r_i} - \bf{cm}$, $\bf{r_i}$ is the coordinate of the $i$th particle, $\bf{cm}$ is the
position of the mass center of the whole group of atoms, $q_i$ is the charge of the $i$th particle,
$\bf{I}$ is the identity matrix and $N$ is the number of atoms.


`multipolemoments`  | Description
------------------- | -------------------------------------
`nstep`             | Interval with which to sample
`indexes`           | Array defining a range of indexes within the molecule 
`index`             | Index of the molecular group
`mol_cm=true`       | Moments w.r.t. the mass center of the whole molecule (instead of the subgroup)

### Electric Multipole Distribution

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
on molecule $a$ and $b$ with charges $q\_i$ in position $\boldsymbol{r}\_i$
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
    u_{\text{dip-dip}}  =  \frac{\boldsymbol{\mu_a}\boldsymbol{\mu_b}   }{ R^3  }
        - \frac{3 (\boldsymbol{\mu_a} \cdot \boldsymbol{R}) ( \boldsymbol{\mu_b}\cdot\boldsymbol{R})  }{R^5} 
$$

$$
    u_{\text{ion-quad}} =  \frac{ q_a \boldsymbol{R}^T \boldsymbol{Q}_b \boldsymbol{R} }{R^5}-\frac{q_a \mbox{tr}(\boldsymbol{Q}_b) }{R^3}+ ...
$$

$$
    u_{\text{total}} =  u_{\text{ion-ion}} + u_{\text{ion-dip}} + u_{\text{dip-dip}} + u_{\text{ion-quad}}
$$

$$
    u_{\text{exact}} =  \sum_i^a\sum_j^b \frac{q_iq_j}{ | \boldsymbol{r_i} - \boldsymbol{r_j} | }
$$


During simulation, the above terms are thermally averaged over angles, co-solute degrees of freedom etc.
Note also that the moments are defined with respect to the _mass_ center, not _charge_ center.
While for globular macromolecules the difference between the two is often small,
the latter is more appropriate and is planned for a future update.

The input keywords are:

`multipoledist`  | Description
---------------- | ----------------------------------------------------
`nstep`          | Interval between samples
`file`           | Output file name
`molecules`      | Array with exactly two molecule names, $a$ and $b$
`dr`             | Distance resolution (Å) along _R_.

### Charge Fluctuations

For a given molecule, this calculates the average charge and standard deviation per atom, 
and the most probable species (atom name) averaged over all present molecules.
A PQR file of a random molecule with average charges and most probable 
atomic species can be saved.

`chargefluctuations` | Description
-------------------- | ----------------------------------------------------------
`nstep`              | Interval between samples
`nskip=0`            | Number of initial steps excluded from the analysis
`molecule`           | Name of molecular group
`pqrfile`            | Output PQR file (optional)
`verbose=True`       | If `True`, add results to general output

### Electric Potential

`electricpotential`   | Description
--------------------- | ------------------------------------------------------
`nstep`               | Interval between samples
`nskip=0`             | Number of initial steps excluded from the analysis
`epsr`                | Dielectric constant
`type`                | Coulomb type, `plain` etc. -- see energies
`structure`           | Either a _filename_ (pqr, aam, gro etc) or a _list_ of positions
`policy=fixed`        | Policy used to augment positions before each sample event, see below
`ncalc`               | Number of potential calculations per sample event
`stride`              | Separation between target points when using `random_walk` or `no_overlap`

This calculates the mean electric potential, $\langle \phi\_i \rangle$ and correlations, $\langle \phi\_1\phi\_2 ...\rangle$
at an arbitrary number of target positions in the simulation cell.
The positions - given via `structure` - can be augmented using a `policy`:

`policy`        | Description
--------------- | ------------------------------------------------------------------
`fixed`         | Expects a list of fixed positions where the potential is measured
`random_walk`   | Assign random position to first; the following targets are randomly placed `stride` distance from previous.
`no_overlap`    | As `random_walk` but with no particle overlap (size defined by `sigma`, see Topology)

Histograms of the correlation and the potentials at the target points are saved to disk.

Example:

~~~ yaml
- electricpotential:
    nstep: 20
    ncalc: 10
    epsr: 80
    type: plain
    policy: random_walk
    stride: 10   # in angstrom
    structure:
      - [0,0,0]  # defines two target points...
      - [0,0,0]  # ...positions are randomly set

## Reaction Coordinate

This saves a given [reaction coordinate](energy.html#reaction-coordinates)
as a function of steps. The generated output `file` has three columns:

1. step number
2. the value of the reaction coordinate
3. the cummulative average of all preceding values.

Optional [gzip compression](https://en.wikipedia.org/wiki/Gzip)
can be enabled by suffixing the filename with `.gz`, thereby reducing the output file size significantly.
The folowing example reports the mass center $z$ coordinate of the first molecule every 100th steps:

~~~ yaml
- reactioncoordinate:
    {nstep: 100, file: cmz.dat.gz, type: molecule, index: 0, property: com_z}
~~~ 

In the next example, the angle between the principal molecular axis and the $xy$-plane
is reported by diagonalising the gyration tensor to find the principal moments:

~~~ yaml
- reactioncoordinate:
    {nstep: 100, file: angle.dat.gz, type: molecule, index: 0, property: angle, dir: [0,0,1]}
~~~ 

### Processing

In the above examples we stored two properties as a function of steps. To join the two
files and calculate the _average angle_ as a function of the mass center coordinate, _z_,
the following python code may be used:

~~~ python
import numpy as np
from scipy.stats import binned_statistic

def joinRC(filename1, filename2, bins):
    x = np.loadtxt(filename1, usecols=[1])
    y = np.loadtxt(filename2, usecols=[1])
    means, edges, bins = binned_statistic(x, y, 'mean', bins)
    return (edges[:-1] + edges[1:]) / 2, means

cmz, angle = joinRC('cmz.dat.gz', 'angle.dat.gz', 100)
np.diff(cmz) # --> cmz resolution; control w. `bins`
~~~

Note that Numpy automatically detects and decompresses `.gz` files.
Further, the command line tools `zcat`, `zless` etc. are useful for handling
compressed files.


## System Sanity

It is wise to always assert that the simulation
is internally sane. This analysis checks the following and aborts if insane:

- all particles are inside the simulation boundaries
- molecular mass centers are correct
- ...more to be added.

To envoke, use for example `- sanity: {nstep: 1}` by default, `nstep=-1`, meaning it will
be run at the end of simulation, only.
This is not a particularly time-consuming analysis and we recommend that it is enabled
for all simulations.


## System Energy

Calculates the energy contributions from all terms in the Hamiltonian and
outputs to a file as a function of steps.
If filename ends with `.csv`, a comma separated value file will be saved,
otherwise a simple space separated file with a single hash commented header line.
To enable GZip compression, suffix the filename with `.gz`.
All units in $k\_BT$.

`systemenergy` | Description
-------------- | -------------------------------------------
`file`         | Output filename (`.dat`, `.csv`, `.dat.gz`)
`nstep`        | Interval between samples


## Perturbations

### Virtual Volume Move

Performs a [virtual volume move](http://doi.org/cppxt6) by
scaling the simulation volume to $V+\Delta V$ along with
molecular mass centers and atomic positions. The excess pressure is evaluated
as a Widom average:

$$
    p^{ex} = \frac{k_BT}{\Delta V} \ln \left\langle e^{-\delta u / k_BT} \right\rangle_{NVT}
$$

If `file` is given, the pressure as a function of steps is written to a (compressed) file.

`virtualvolume`     | Description
------------------- | -------------------------------------
`dV`                | Volume perturbation (Å³)
`nstep`             | Interval between samples
`file`              | Optional output filename (`.dat`, `.dat.gz`)
`scaling=isotropic` | Volume scaling method (`isotropic`, `xy`, `z`)

By default, the volume is isotropically scaled, but for more advanced applications of
volume perturbations - pressure tensors, surface tension etc., see [here](http://doi.org/ckfh).
If a non-isotropic scaling is used, an extra column will be added to the output
`file` containing the change in area (`xy`) or length (`z`).
See also the documentation for the Monte Carlo _Volume move_.


### Virtual Translate Move

Performs a virtual displacement, $dL$, of a single `molecule` in
the direction `dir` and measure the force by perturbation,

$$
    f = \frac {k_BT \ln \langle \exp{\left (-dU/k_BT \right )} \rangle_0 }{ dL }
$$

`virtualtranslate` | Description
------------------ | ---------------------------------------------------------------
`molecule`         | Molecule name; only _one_ of these is allowed in the system
`dL`               | Displacement (Å)
`dir=[0,0,1]`      | Displacement direction (length ignored)
`nstep`            | Interval between samples
`file`             | Optional output filename for writing data as a function of steps (`.dat|.dat.gz`)


### Widom Insertion

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
inserted group. This is useful for example to calculate the excess activity
coefficient of a neutral salt pair.
Upon insertion, random positions and orientations are
generated. For use with rod-like particles on surfaces, the `absz`
keyword may be used to ensure orientations on only one
half-sphere.

Exactly _one inactive_ `molecule` must be added to the simulation using the `inactive`
keyword when inserting the initial molecules in the [topology](topology).

`widom`       | Description
------------- | -----------------------------------------
`molecule`    | Name of _inactive_ molecule to insert (atomic or molecular)
`ninsert`     | Number of insertions per sample event
`dir=[1,1,1]` | Inserting directions
`absz=false`  | Apply `std::fabs` on all z-coordinates of inserted molecule
`nstep`       |  Interval between samples

## Positions and Trajectories

### Save State

`savestate`        |  Description
------------------ | ------------------------------------------------------------------------------------------
`file`             |  File to save; format detected by file extension: `pqr`, `aam`, `gro`, `xyz`, `json`/`ubj`
`saverandom=false` |  Save the state of the random number generator
`nstep=-1`         |  Interval between samples; if -1 save at end of simulation
`convert_hexagon`  |  Convert hexagonal prism to space-filling cuboid; `pqr` only (default: false)

Saves the current configuration or the system state to file. For grand canonical
simulations, the PQR file format sets charges and radii of inactive particles to zero
and positions them in one corner of the box.

If the suffix is `json` or `ubj` ([binary](http://ubjson.org)), a single 
state file that can be used to restart the simulation is saved
with the following information:

- topology: atom, molecule, and reaction definitions
- particle and group properties incl. positions
- geometry
- state of random number generator (if `saverandom=true`)

If `nstep` is greater than zero, the output filename will be tagged
with the current step count.


### Space Trajectory (experimental)

Save all particle and group information to a compressed, binary trajectory format.
The following properties are saved:

 - all particle properties (id, position, charge, dipole etc.)
 - all group properties (id, size, capacity etc.)
 - todo: geometry, energy

The file suffix must be either `.traj` (uncomressed) or `.ztraj` (compressed).
For the latter, the file size is reduced by roughly a factor of two using zlib
compression.

`spacetraj`  | Description
------------ | ---------------------------------------
`file`       | Filename of output .traj/.ztraj file
`nstep`      | Interval between samples.


### XTC trajectory

Generates a Gromacs XTC trajectory file with particle positions and box
dimensions as a function of steps. Both _active_ and _inactive_ atoms are
saved.

`xtcfile`      |  Description
-------------- | ---------------------------------------------------------
`file`         |  Filename of output xtc file
`nstep`        |  Interval between samples.
`molecules=*`  |  Array of molecules to save (default: all)


### Charge-Radius trajectory

Most trajectory file formats do not support a fluctuating number
of particles. For each `nstep`, this analysis files charge and
radius information for all particles.
Inactive particles are included with _zero_ charge and radius.

Using a helper script for VMD (see `scripts/`) this information
can be loaded to visualise flutuating charges and or number of particles.
The script should be sourced from the VMD console after loading the trajectory,
or invoked when launching VMD:

~~~ bash
vmd confout.pqr traj.xtc -e scripts/vmd-qrtraj.tcl
~~~

`qrfile`          |  Description
----------------- | -----------------------------------
`file=qrtraj.dat` |  Output filename (.dat, .gz)
`nstep`           |  Interval between samples
