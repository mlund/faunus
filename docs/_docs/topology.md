<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
});
</script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

# Topology

The topology describes atomic and molecular properties as well as processes and reactions.

## Global Properties

The following keywords control temperature, simulation box size etc., and must be
placed outer-most in the input file.

~~~ yaml
temperature: 298.15  # system temperature (K)
geometry:
  type: cuboid       # Cuboidal simulation container
  length: [40,40,40] # cuboid dimensions (array or number)
mcloop:              # number of MC steps (macro × micro)
  macro: 5           # Number of outer MC steps
  micro: 100         # Number of inner MC steps; total = 5 × 100 = 500
random:              # seed for pseudo random number generator
  seed: fixed        # "fixed" (default) or "hardware" (non-deterministic)
~~~

### Geometry

Below is a list of possible geometries, specified by `type`, for the simulation container,
indicating if and in which directions periodic boundary conditions (PBC) are applied.
Origin ($0,0,0$) is always placed in the geometric _center_ of the simulation container.
Particles are always kept inside the simulation container with an external
potential that is zero if inside; infinity if outside.

`geometry` | PBC      | Required keywords
---------- | -------- | --------------------------------------
`type`:    |          |
`cuboid`   | $x,y,z$  | `length` (array or single float)
`slit`     | $x,y$    | `length` (array or single float)
`hexagonal`| $x,y$    | `radius` (inscribed/inner), `length` (along _z_)
`cylinder` | $z$      | `radius`, `length` (along _z_)
`sphere`   | none     | `radius`

### Simulation Steps

The variables `macro` and `micro` are positive integers and their product
defines the total number simulations steps.
In each step a random Monte Carlo move is drawn from a weighted distribution.
For each `macro` step, all analysis methods are, if befitting, instructed to
flush buffered data to disk and may also trigger terminal output.
For this reason `macro` is typically set lower than `micro`.

## Atom Properties

Atoms are the smallest possible particle entities with properties defined below.

`atomlist`    | Description
------------- | ------------------------------------------------------
`activity=0`  | Chemical activity for grand canonical MC [mol/l]
`pactivity`   | −log10 of chemical activity (will be converted to activity)
`alphax=0`    | Excess polarizability (unit-less)
`dp=0`        | Translational displacement parameter [Å]
`dprot=0`     | Rotational displacement parameter [radians]
`eps=0`       | Lennard-Jones/WCA energy parameter [kJ/mol]
`mu=[0,0,0]`  | Dipole moment vector [eÅ]
`mulen=|mu|`  | Dipole moment scalar [eÅ]
`mw=1`        | Molecular weight [g/mol]
`q=0`         | Valency or partial charge number [$e$]
`r=0`         | Radius = `sigma/2` [Å]
`sigma=0`     | `2r` [Å] (overrides radius)
`tension=0`   | Surface tension [kJ/mol/Å$^2$]
`tfe=0`       | Transfer free energy [kJ/mol/Å$^2$/M]

A filename (`.json`) may be given instead of an atom definition to load
from an external atom list. Atoms are loaded in the given order, and if it occurs
more than once, the latest entry is used.

Example:

~~~ yaml
atomlist:
  - Na: {q:  1.0, sigma: 4, eps: 0.05, dp: 0.4}
  - Ow: {q: -0.8476, eps: 0.65, sigma: 3.165, mw: 16}
  - my-external-atomlist.json
  - ...
~~~


## Molecule Properties

A molecule is a collection of atoms, but need not be associated
as real molecules. Two particular modes can be specified:

1. If `atomic=true` the atoms in the molecule are unassociated and is
   typically used to define salt particles or other non-aggregated
   species. No structure is required, and the molecular center of mass (COM) is
   unspecified.

2. If `atomic=false` the molecule resembles a real molecule and a structure or trajectory
   is _required_.

Properties of molecules and their default values:

`moleculelist`          | Description
------------------      | -----------------------------------------------------
`activity=0`            | Chemical activity for grand canonical MC [mol/l]
`atomic=false`          | True if collection of atomic species, salt etc.
`atoms=[]`              | Array of atom names; required if `atomic=true`
`bondlist`              | List of _internal_ bonds (harmonic, dihedrals etc.)
`compressible=false`    | If true, molecular internal coordinates are scaled upon volume moves
`ensphere=false`        | Radial rescale of positions to sphere w. radius of average radial distance from COM (stored in 1st atom which is a dummy)
`excluded_neighbours=0` | Generate an `exclusionlist` from the bonded interaction: Add all atom pairs which are `excluded_neighbours` or less bonds apart
`exclusionlist`         | List of _internal_ atom pairs which nonbonded interactions are excluded
`implicit=false`        | Mark as implicit for reactive Monte Carlo schemes
`insdir=[1,1,1]`        | Insert directions are scaled by this
`insoffset=[0,0,0]`     | Shifts mass center after insertion
`keeppos=false`         | Keep original positions of `structure`
`keepcharges=true`      | Keep charges of `structure` (aam/pqr files) and `traj` even if mismatch with those in `atomlist`
`rigid=false`           | Set to true for rigid molecules. Affects energy evaluation.
`rotate=true`           | If false, the original structure will not be rotated upon insertion
`structure`             | Structure file or direct information; required if `atomic=false`
`to_disk=false`         | Save initial structure to `{name}-initial.pqr`; for molecular groups only
`traj`                  | Read conformations from PQR trajectory. Cannot be used w. `structure`; see also `keepcharges`
`trajweight`            | One-column file with relative weights for each conformation. Must match frames in `traj` file.
`trajcenter=false`      | Move CM of conformations to the origin assuming whole molecules

Example:

~~~ yaml
moleculelist:
  - salt: {atoms: [Na,Cl], atomic: true}
  - water:
      structure: water.xyz
      bondlist:
        - harmonic: {index: [0,1], k: 100, req: 1.5}
        - ...
  - carbon_dioxide:
      structure:
        - O: [-1.162,0,0]
        - C: [0,0,0]
        - O: [1.162,0,0]
      bondlist:
        - harmonic: {index: [1,0], k: 8443, req: 1.162}
        - harmonic: {index: [1,2], k: 8443, req: 1.162}
        - harmonic_torsion: {index: [0,1,2], k: 451.9, aeq: 180}
      excluded_neighbours: 2 # generates an exclusionlist as shown below
      exclusionlist: [ [0,1], [1,2], [0,2] ] # redundant in this topology
  - ...
~~~

### Structure Loading Policies

When giving structures using the `structure` keyword, the following policies apply:

- `structure` can be a file name: `file.@` where `@=xyz|pqr|aam|gro`
- `structure` can be an _array_ of atom names and their positions:
  `- Mg: [2.0,0.1,2.0]`
- `structure` can be a [FASTA sequence](https://en.wikipedia.org/wiki/FASTA_format):
  `{fasta: AAAAAAAK, k: 2.0; req: 7.0}` which generates
  a linear chain of harmonically connected atoms.
  FASTA letters are translated into three letter residue names which _must_ be defined
  in `atomlist`.
  Special letters: `n=NTR`, `c=CTR`, `a=ANK`.
  Instead of a sequence, `fasta` may be a _filename_ from which the first
  sequence is extracted. The filename must end with `.fasta`.
- Radii in files are _ignored_; `atomlist` definitions are used.
  A warning is issued if radii/charges differ in files and `atomlist`.
- By default, charges in files are _used_; `atomlist` definitions are ignored.
  Use `keepcharges=False` to override.
- If the structure file contains box size information, this will be _ignored_.

### Nonbonded Interaction Exclusion

Some nonbonded interactions between atoms within a molecule may be excluded in the topology.
Force fields almost always exclude nonbonded interactions between directly bonded atoms. However
other nonbonded interactions may be excluded as well; refer to your force field parametrization.
If a molecule contains overlapping hard spheres, e.g., if the bond length is shorter than
the spheresʼ diameter, it is necessary to exclude corresponding nonbonded interactions to avoid
infinite energies.

The excluded nonbonded interactions can be given as an explicit list of atom pairs `excludelist`,
or they can be deduced from the moleculeʼs topology using the `excluded_neighbours=n` option:
If the atoms are `n` or less bonds apart from each other in the molecule, the nonbonded interactions
between them are excluded. Both options `excluded_neighbours` and `exclusionlist` can be used together
making a union.

## Initial Configuration

Upon starting a simulation, an initial configuration is required and must be
specified in the section `insertmolecules` as a list of valid molecule names.
Molecules are inserted in the given order and may be `inactive`, meaning that
they are not present in the simulation cell, but available as a reservoir for
e.g. grand canonical moves.
If a group is marked `atomic`, its `atoms` are inserted `N` times.

Example:

~~~ yaml
insertmolecules:
  - salt:  { molarity: 0.1 }
  - water: { N: 256, inactive: 2 }
~~~

The following keywords for each molecule type are available:

`insertmolecules`    | Description
-------------------- | ---------------------------------------
`N`                  | Number of molecules to insert
`molarity`           | Insert molecules to reach molarity
`inactive`           | Number of inserted molecules to deactivate; set to `true` for all
`positions`          | Load positions from file (`aam`, `pqr`, `xyz`)
`translate=[0,0,0]`  | Displace loaded `positions` with vector

A filename with positions for the `N` molecules can be given with `positions`.
The file must contain exactly `N`-times molecular
positions that must all fit within the simulation box. Only _positions_ from
the file are copied; all other information is ignored.

For `implicit` molecules, only `N` should be given and the molecules are never
inserted into the simulation box.

The `molarity` keyword is an alternative to `N` and uses the initial
volume to calculate the number of molecules to insert. `N` and
`molarity` are mutually exclusive.


### Overlap Check

Random insertion is repeated until there is no overlap with the simulation
container boundaries. Overlap between particles is ignored and for
i.e. hard-sphere potentials the initial energy may be infinite.


## Equilibrium Reactions

Faunus supports density fluctuations, coupled to chemical equilibria with
explicit and/or implicit particles via their chemical potentials as
defined in the `reactionlist` detailed below, as well as in `atomlist` and
`moleculelist`. The level of flexibility is high and reactions can be
freely composed.

The move involves deletion and insertion of reactants and products and it is
therefore important that simulations are started with a sufficiently high number of
initial molecules in `insertmolecules`.
If not, the `rcmc` move will attempt to issue warnings with suggestions how to fix it.

### Reaction format

The initial key describes a transformation of reactants (left of `=`)
into products (right of `=`) that may be a mix of atomic and molecular species.

- all species, `+`, and `=` must be surrounded by white-space
- atom and molecule names cannot overlap
- species can be repeated to match the desired stoichiometry, _e.g._ `A + A = C`

Available keywords:

`reactionlist`  | Description
--------------- | ---------------------------------------------------------------
`lnK`/`pK`      | Molar equilibrium constant either as $\ln K$ or $-\log\_{10}(K)$
`neutral=false` | If true, only neutral molecules participate in the reaction

The `neutral` keyword is needed for molecular groups containing titratable atoms. If `neutral` is set to true,
the activity of the neutral molecule should be specified in `moleculelist`.


### Example: Grand Canonical Salt Particles

This illustrates how to maintain constant chemical potential of salt ions:

~~~ yaml
atomlist:
  - na: {q: 1.0, ...}  # note that atom names must differ
  - cl: {q: -1.0, ...} # from molecule names
moleculelist:
  - Na+: {atoms: [na], atomic: true, activity: 0.1}
  - Cl-: {atoms: [cl], atomic: true, activity: 0.1}
reactionlist:
  - "= Na+ + Cl+": {} # note: molecules, not atoms
moves:
  - rcmc: {} # activate speciation move
~~~

The same setup can be used also for molecular molecules, _i.e._ molecules with `atomic: false`.


### Example: Acid-base titration with _implicit_ protons

An _implicit_ atomic reactant or product is included in the reaction but not 
explicitly in the simulation cell.
Common use-cases are acid-base equilibria where the proton concentration is often very low:

~~~ yaml
atomlist:
  - H+: {implicit: true, activity: 0.00001} # pH 5
  - COO-: {q: -1.0, ...}
  - COOH: {q: 0.0, ...}
reactionlist:
  - "COOH = COO- + H+": {pK: 4.8} # not electroneutral!
~~~

where we set `pK` equal to the `pKa`:
$$
K\_a = \frac{ a_{\mathrm{COO^-}} a_{\mathrm{H^+}} }{ a_{\mathrm{COOH}} }.
$$
To simulate at a given constant pH, H+ is specified as an implicit atom of activity $10^{-\mathrm{pH}}$ and the equilibrium 
is modified accordingly (in this case $K$ is divided by $a_{\mathrm{H^+}}$). 
It is important to note that this reaction violates _electroneutrality_ and should be used
only with Hamiltonians where this is allowed. This could for example be in systems with salt screened
Yukawa interactions. 


### Example: Acid-base titration coupled with Grand Canonical Salt

To respect electroneutrality when swapping species, we can associate the titration move with
an artificial insertion or deletion of salt ions. These ions should be present under constant
chemical potential and we therefore couple to a grand canonical salt bath:

~~~ yaml
atomlist:
  - H+: {implicit: true, activity: 0.00001} # pH 5
  - COO-: {q: -1.0, ...}
  - COOH: {q: 0.0, ...}
  - na: {q: 1.0, ...}
  - cl: {q: -1.0, ...}
moleculelist:
  - Na+: {atoms: [na], atomic: true, activity: 0.1}
  - Cl-: {atoms: [cl], atomic: true, activity: 0.1}
reactionlist:
  - "COOH + Cl- = COO- + H+": {pK: 4.8} # electroneutral!
  - "COOH = Na+ + COO- + H+": {pK: 4.8} # electroneutral!
  - "= Na+ + Cl-": {} # grand canonical salt
~~~

For the first reaction, $K$ is divided by both $a_{\mathrm{H^+}}$ and $a_{\mathrm{Cl^-}}$, so that the final equilibrium constant
used by the speciation move is
$$
K' = \frac{K\_a}{a_{ \mathrm{H^+} } a_{ \mathrm{Cl^-} } } = \frac{ a_{\mathrm{COO^-}} }{ a_{\mathrm{COOH}} a_{\mathrm{Cl^-}} }.
$$
In an ideal system, the involvement of Na or Cl in the acid-base reaction is inconsequential for the equilibrium,
since the Grand Canonical ensemble ensures constant salt activity.


### Example: Precipitation of Calcium Hydroxide using _implicit_ molecules

Here we introduce a solid phase of Ca(OH)₂ and its solubility product
to predict the amount of dissolved calcium and hydroxide ions. Note that
we start from an empty simulation box (both ions are inactive) and the solid
phase is treated _implicitly_, _i.e._ it never inters the simulation box.
If a forward reaction is made, one implicit molecule (out 200 in total) will
be converted into explicit molecules.
Once all 200 have been consumed, only backward reactions are possible.
Additional, coupled reactions can be introduced to study complex
equilibrium systems under influence of intermolecular interactions.

~~~ yaml
moleculelist:
    - Ca(OH)2: {implicit: true} # this molecule is implicit
    - Ca++: {atoms: [ca++], atomic: true}
    - OH-: {atoms: [oh-], atomic: true}
insertmolecules:
    - Ca++: {N: 200, inactive: true}
    - OH-: {N: 400, inactive: true}
    - Ca(OH)2: {N: 200} # not actually inserted!
reactionlist:
    - "Ca(OH)2 = Ca++ + OH- + OH-": {pK: 5.19}
~~~


### Example: Swapping between molecular conformations

The following can be used to alternate between different molecular conformations.
When swapping, the mass center position and orientation are randomly generated.

~~~ yaml
moleculelist:
  - A: {structure: "conformationA.xyz", ...}
  - B: {structure: "conformationB.xyz", ...}
reactionlist:
  - "A = B": {lnK: 0.69} # K=2, "B" twice as likely as "A"
~~~

### Internal degrees of freedom (experimental)

If a fluctuating molecule has internal degreees of freedom, the internal bond energy is included
as a bias so that the internal state does not affect the acceptance.
To disable this behavior, a minor code modification is currently required (see `MolecularGroupDeActivator::apply_bond_bias`).

