---
---
<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
});
</script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>
[![Edit](https://img.shields.io/badge/Github-Improve_this_page-orange.svg)]({{site.github.repository_url}}/blob/master/docs/{{page.path}})

# Topology

The topology describes atomic and molecular properties as well as processes and reactions.

## Global Properties

The following keywords control temperature, simulation box size etc., and must be
placed outer-most in the input file.

### System Temperature

-------------- | -----------------------
`temperature`  | Temperature [K]

### Simulation Container

If a single number is given, a cube is generated, while an array
of size three is used to set different side-lengths.

`geometry`     | Description
-------------- | -----------------------------------------------------
`length`       | Cuboidal side-length(s) [Å] (float or array with $L_x, L_y, L_z$)

### Number of Monte Carlo Loops

The total number of MC steps is `macro` x `micro`.

`mcloop`       | Description
-------------- | ---------------------
`macro`        | Number of macro loops (integer)
`micro`        | Number of micro loops (integer)

### Random Number Generator

By default a deterministic sequence is generated, while
the `hardware` seed attempts to use a hardware seed to provide
a _non-deterministric_ sequence.

`random`       | Description 
-------------- | -----------------------------
`seed`         | `fixed` (default) or `hardware`

## Atom Properties

Atoms are the smallest possible particle entities with properties defined below.

`atomlist`    | Description
------------- | ------------------------------------------------------
`activity=0`  | Chemical activity for grand canonical MC [mol/l]
`alphax=0`    | Excess polarizability (unit-less)
`dp=0`        | Translational displacement parameter [Å]
`dprot=0`     | Rotational displacement parameter [degrees] (will be converted to radians)
`eps=0`       | Epsilon energy scaling commonly used for Lennard-Jones interactions etc. [kJ/mol]
`mu=[0,0,0]`  | Dipole moment vector [Debye]
`mw=1`        | Molecular weight [g/mol]
`q=0`         | Valency / partial charge number [$e$]
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
   typically used to defined salt particles or any other non-aggregated
   species. No structure is required, and the molecular center of mass (COM) is
   unspecified.

2. If `atomic=false` the molecule resembles a real molecule and a structure
   is _required_.

Properties of molecules and their default values:

`moleculelist`      | Description
------------------- | -------------------------------------------------
`activity=0`        | Chemical activity for grand canonical MC [mol/l]
`atomic=false`      | True if collection of atomic species, salt etc.
`atoms=[]`          | Array of atom names - required if `atomic=true`
`implicit=false`    | If this species is implicit in GCMC schemes
`insdir=[1,1,1]`    | Insert directions are scaled by this
`insoffset=[0,0,0]` | Shifts mass center after insertion
`keeppos=false`     | Keep original positions of `structure`
`structure`         | Structure file (`.pqr, .aam, .xyz`) - required if `atomic=false`
`bondlist`          | List of _internal_ bonds (harmonic, dihedrals etc.)
`traj`              | Read conformations from PQR trajectory (`structure` will be ignored)
`trajweight`        | One column file w. relative weights for each conformation. Must match frames in `traj` file.
`trajcenter`        | CM of conformations to origo assuming whole molecules (default: `false`)

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

If the file given by `structure` contains charges or radii (i.e. pqr, aam files), this information
is used over data from `atomlist`.


### Initial Configuration

Upon starting a simulation, an initial configuration is required and must be
specified in the section `insertmolecules` as a list of valid molecule names.
Molecules are inserted in the given order and may be `inactive`
which is used by some analysis and ensembles.
If a group is marked `atomic`, its `atoms` will be inserted `N` times.

Example:

~~~ yaml
insertmolecules:
    - salt:  { N: 10 }
    - water: { N: 256 }
    - water: { N: 1, inactive: true }
~~~

The following keywords are available:

Keyword             | Description
------------------- | -------------------------------------------------
`N`                 | Number of molecules to insert
`inactive=false`    | Deactivates inserted molecules
`positions`         | Load positions from file (`aam`, `pqr`, `xyz`)
`translate=[0,0,0]` | Displace loaded `positions` with vector

A filename with positions for the N molecules can be given with `positions`.
The file must contain exactly N-times molecular
positions that must all fit within the simulation box. Only positions from
the file are copied while all other information is ignored.

## Equilibrium Reactions

Faunus supports density fluctuations, coupled to chemical equilibria with
explicit and/or implicit particles via their chemical potentials as
defined in the `reactionlist` detailed below, as well as in `atomlist` and
`moleculelist`.

~~~ yaml
reactionlist:
    - "AH = A + H": { pK: 4.8, canonic: true }
    - "Mg(OH)2 = Mg + OH + OH": { lnK: -25.9 }
~~~

The initial string describes a transformation of reactants (left of `=`)
into products (right of `=`) that may be a mix of atomic and molecular species.

Note that:

- all species, `+`, and `=` must be surrounded by white-space
- atom and molecule names cannot overlap
- you may repeat species to match desired stoichiometry

Available keywords:

`reactionlist`  | Description
--------------- | ---------------------------------------------------------------
`lnK`/`pK`      | Molar equilibrium constant either as $\ln K$ or $-\log_{10}(K)$

