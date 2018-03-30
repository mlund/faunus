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

### Initial Configuration

Upon starting a simulation, an initial configuration is required and must be
specified in the section `insertmolecules` as a list of valid molecule names.
Molecules are inserted in the given order and may be `inactive`
which is used by some analysis and ensembles.
If a group is marked `atomic`, its `atoms` will be inserted `N` times.

Example:

~~~ yaml
insertmolecules:
    - salt: { N: 10 }
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

## Processes

`processlist`   | Description
--------------- | ----------------------------------------------
`process`       | Process involving molecular groups (string)
`K`/`pK`        | Molar equilibrium constant or minus log thereof

The `process` string describes a transformation of molecular reactants (left of =)
into molecular products (right of =). A trailing tilde (~) denotes that
the species is to be treated _implicitly_.
In the example below we assume that `HA`, `H`, and `A` have
been defined in `moleculelist`:

~~~ yaml
processlist:
    - { process: "HA = H~ + A", pK=4.8 }
~~~
