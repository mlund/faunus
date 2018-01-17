---
---
<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
});
</script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

# Topology

The topology describes atomic and molecular properties and is defines at the top-level input:

~~~ yaml
atomlist:
    - Na: {q:  1.0, sigma: 4, eps: 0.05, dp: 0.4}
    - Ow: {q: -0.8476, eps: 0.65, sigma: 3.165, mw: 16}
    - ...
moleculelist:
    - salt: {Ninit: 10, atoms: [Na,Cl], atomic: true}
    - water:
        Ninit: 256
        structure: water.xyz
        bondlist:
            - harmonic: {index: [0,1], k=100; req=1.5}
            - ...
    - ...
~~~

## Atom Properties

Atoms are the smallest possible particle entities and their properties are defined in the
Table below.

`atomlist`    | Description
------------- | ------------------------------------------------------
`activity=0`  | Chemical activity for grand canonical MC [mol/l]
`dp=0`        | Translational displacement parameter [Å]
`dprot=0`     | Rotational displacement parameter [degrees] (will be converted to radians)
`eps=0`       | Epsilon energy scaling commonly used for Lennard-Jones interactions etc. [kJ/mol]
`mu=[0,0,0]`  | Dipole moment vector [Debye]
`Ninit=0`     | Initial number of atoms (used by `MoleculeData` to insert atoms
`mw=1`        | Molecular weight [g/mol]
`q=0`         | Valency / partial charge number [$e$]
`r=0`         | Radius = `sigma/2` [Å]
`sigma=0`     | `2r` [Å] (overrides radius)
`tension=0`   | Surface tension [kJ/mol/Å$^2$]
`tfe=0`       | Transfer free energy [kJ/mol/Å$^2$/M]

## Molecule Properties

A molecule is a collection of atoms, but they need not be associated
as real molecules. Two particular modes can be specified.

1. If `atomic=true` the atoms in the molecule are unassociated and is
   typically used to defined salt particles or any other non-aggregated
   species. No structure is required, and the molecular center of mass (COM) is
   unspecified. `Ninit` is used to insert _N_-times the number of
   atoms defined in `atoms`.

2. If `atomic=false` the molecule resembles a real molecule and a structure
   is _required_. `Ninit` is used to insert _N_ molecules.

Properties of molecules and their default values:

`moleculelist`      | Description
------------------- | -------------------------------------------------
`activity=0`        | Chemical activity for grand canonical MC [mol/l]
`atomic=false`      | True if collection of atomic species, salt etc.
`atoms=[]`          | Array of atom names - required if `atomic==true`
`implicit=false`    | If this species is implicit in GCMC schemes
`insdir=[1,1,1]`    | Insert directions are scaled by this
`insoffset=[0,0,0]` | Shifts mass center after insertion
`keeppos=false`     | Keep original positions of `structure`
`Ninactive=0`       | Deactivates `Ninactive` of the inserted molecules
`Ninit=0`           | How many molecules to insert
`structure`         | Structure file (`.pqr|.aam|.xyz`) - required if `atomic=false`
`bondlist`          | List of _internal_ bonds (harmonic, dihedrals etc.)

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
