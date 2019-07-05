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
mcloop:              # number of MC steps (macro x micro)
  macro: 5           # Number of outer MC steps
  micro: 100         # Number of inner MC steps; total = 5 x 100 = 5000
random:              # seed for pseudo random number generator
  seed: fixed        # "fixed" (default) or "hardware" (non-deterministic)
~~~

### Geometry

Below is a list of possible geometries, specified by `type`, for the simulation container,
indicating if and in which directions periodic boundary conditions (PBC) are applied.
Origin ($0,0,0$) is always placed in the geometric _center_ of the simulation container.

`geometry` | PBC      | Required keywords
---------- | -------- | --------------------------------------
`type`:    |          |
`cuboid`   | $x,y,z$  | `length` (array or single float)
`slit`     | $x,y$    | `length` (array or single float)
`hexagonal`| $x,y$    | `radius` (inscribed/inner), `length` (along _z_)
`cylinder` | $z$      | `radius`, `length` (along _z_)
`sphere`   | none     | `radius`

## Atom Properties

Atoms are the smallest possible particle entities with properties defined below.

`atomlist`    | Description
------------- | ------------------------------------------------------
`activity=0`  | Chemical activity for grand canonical MC [mol/l]
`pactivity`   | -log10 of chemical activity (will be converted to activity)
`alphax=0`    | Excess polarizability (unit-less)
`dp=0`        | Translational displacement parameter [Å]
`dprot=0`     | Rotational displacement parameter [degrees] (will be converted to radians)
`eps=0`       | Epsilon energy scaling commonly used for Lennard-Jones interactions etc. [kJ/mol]
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

`moleculelist`      | Description
------------------- | -------------------------------------------------
`activity=0`        | Chemical activity for grand canonical MC [mol/l]
`atomic=false`      | True if collection of atomic species, salt etc.
`atoms=[]`          | Array of atom names - required if `atomic=true`
`bondlist`          | List of _internal_ bonds (harmonic, dihedrals etc.)
`implicit=false`    | If this species is implicit in GCMC schemes
`insdir=[1,1,1]`    | Insert directions are scaled by this
`insoffset=[0,0,0]` | Shifts mass center after insertion
`keeppos=false`     | Keep original positions of `structure`
`keepcharges=true`  | Keep original charges of `structure` (aam/pqr files)
`rigid=false`       | Set to true for rigid molecules. Affects energy evaluation.
`rotate=true`       | If false, the original structure will not be rotated upon insertion
`structure`         | Structure file or direct information - required if `atomic=false`
`traj`              | Read conformations from PQR trajectory (`structure` will be ignored)
`trajweight`        | One column file w. relative weights for each conformation. Must match frames in `traj` file.
`trajcenter=false`  | Move CM of conformations to origo assuming whole molecules (default: `false`)

Example:

~~~ yaml
moleculelist:
  - salt: {atoms: [Na,Cl], atomic: true}
  - water:
      structure: water.xyz
      bondlist:
        - harmonic: {index: [0,1], k: 100, req: 1.5}
        - ...
  - ...
~~~

### Structure Loading Policies

When giving structures using the `structure` keyword, the following policies apply:

- `structure` can be a file name: `file.@` where `@=xyz|pqr|aam`
- `structure` can be an _array_ of atom names and their positions:
  `- Mg: [2.0,0.1,2.0]`
- `structure` can be a [FASTA sequence](https://en.wikipedia.org/wiki/FASTA_format):
  `{fasta: [AAAAAAAK], k: 2.0; req: 7.0}` which generates
  a linear chain of harmonically connected atoms.
  FASTA letters are translated into three letter residue names which _must_ be defined
  in `atomlist`.
  Special letters: `n=NTR`, `c=CTR`, `a=ANK`.
- Radii in files are _ignored_; `atomlist` definitions are used.
- By default, charges in files are _used_; `atomlist` definitions are ignored.
  Use `keepcharges=False` to override.
- A warning is issued if radii/charges differ in files and `atomlist`.
- Box dimensions in files are ignored.

## Initial Configuration

Upon starting a simulation, an initial configuration is required and must be
specified in the section `insertmolecules` as a list of valid molecule names.
Molecules are inserted in the given order and may be `inactive`.
If a group is marked `atomic`, its `atoms` is inserted `N` times.

Example:

~~~ yaml
insertmolecules:
  - salt:  { N: 10 }
  - water: { N: 256 }
  - water: { N: 1, inactive: true }
~~~

The following keywords for each molecule type are available:

`insertmolecules`    | Description
-------------------- | ---------------------------------------
`N`                  | Number of molecules to insert
`inactive=false`     | Deactivates inserted molecules
`positions`          | Load positions from file (`aam`, `pqr`, `xyz`)
`translate=[0,0,0]`  | Displace loaded `positions` with vector

A filename with positions for the `N` molecules can be given with `positions`.
The file must contain exactly N-times molecular
positions that must all fit within the simulation box. Only _positions_ from
the file are copied; all other information is ignored.

### Overlap Check

Random insertion is repeated until there is no overlap with the simulation
container boundaries. Overlap between particles is ignored and for
i.e. hard-sphere potentials the initial energy may be infinite.


## Equilibrium Reactions

Faunus supports density fluctuations, coupled to chemical equilibria with
explicit and/or implicit particles via their chemical potentials as
defined in the `reactionlist` detailed below, as well as in `atomlist` and
`moleculelist`.
The initial key describes a transformation of reactants (left of `=`)
into products (right of `=`) that may be a mix of atomic and molecular species.

An implicit reactant or product is an atom which is included in the equilibrium constant but it is not represented
explicitly in the simulation cell.
A common example is the acid-base equilibrium of the aspartic acid (treated here as atomic particle):

~~~ yaml
reactionlist:
  - "HASP = ASP + H": { pK: 4.0 }
~~~

where H is defined as _implicit_ in the `atomlist`:

~~~ yaml
  - H: { implicit: true, activity: 1e-7 }
~~~

and we set `pK` equal to the `pKa`, i.e.,
$$
K_a = \frac{ a_{\mathrm{ASP}} a_{\mathrm{H}} }{ a_{\mathrm{HASP}} }.
$$
To simulate at a given constant pH, H is specified as an implicit atom of activity $10^{-\mathrm{pH}}$ and the equilibrium 
is modified accordingly (in this case K is divided by $a_{\mathrm{H}}$). 
An acid-base equilibrium, or any other single-atom ID transformation (see the Move section), can also be coupled with the insertion/deletion
of a molecule. For example, 

~~~ yaml
reactionlist:
  - "HASP + Cl = ASP + H": { pK: 4.0 }
  - "= Na + Cl": { }
~~~

where Na and Cl are included in the `moleculelist` as

~~~ yaml
  - Cl: {atoms: [cl], atomic: true, activity: 0.1 } 
  - Na: {atoms: [na], atomic: true, activity: 0.1 } 
~~~

In this case K is both divided by $a_{\mathrm{H}}$ and $a_{\mathrm{Cl}}$, so that the actual equilibrium constant used by the speciation move is

$$
K' = \frac{K_a}{a_{ \mathrm{H} } a_{ \mathrm{Cl} } } = \frac{ a_{\mathrm{ASP}} }{ a_{\mathrm{HASP}} a_{\mathrm{Cl}} }.
$$

In an ideal system, the involvement of Cl in the acid-base reaction does not affect the equilibrium since the grand canonical ensemble
ensures that the activity of Cl matches its concentration in the simulation cell.

Reaction format:

- all species, `+`, and `=` must be surrounded by white-space
- atom and molecule names cannot overlap
- you may repeat species to match the desired stoichiometry

Available keywords:

`reactionlist`  | Description
--------------- | ---------------------------------------------------------------
`lnK`/`pK`      | Molar equilibrium constant either as $\ln K$ or $-\log_{10}(K)$

