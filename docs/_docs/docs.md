---
permalink: /docs/
---
<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
});
</script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>
[![Edit](https://img.shields.io/badge/Github-Improve_this_page-orange.svg)]({{site.github.repository_url}}/blob/master/docs/{{page.path}})

# Introduction

## About

Faunus is a general Monte Carlo simulation code, designed to be flexible, easy
to use and to modify. The code is written in C++14 and (limited) Python bindings
are available.
The development is a **team effort** with, in chronological order,
many valiant contributions from:

_Mikael Lund, Bjorn Persson, Martin Trulsson,
Ondrej Marsalek, Christophe Labbez, Andre Teixeira,
Anil Kurut, Chris Evers, Magnus Ullner,
Robert Vacha, Axel Thuresson, Bjorn Stenqvist,
Joao Henriques, Alexei Abrikossov, Giulio Tesei,
Lukas Sukenik, Coralie Pasquier, Niels Kouwenhoven_

### Supporting Us

Should you find Faunus useful, do consider supporting us by crediting:

- Stenqvist _et al._ [_Molecular Simulation 2013, 39:1233_](http://dx.doi.org/10/nvn)
- Lund, M., Persson, B., Trulsson, M. [_Source Code Biol. Med., 2008, 3:1_](http://dx.doi.org/10/dfqgch)

## Quick Start

Below is an example, `input.yml`, for a Metropolis Monte Carlo simulation
of charged Lennard-Jones particles in a cubic PBC box. Running it with

~~~ bash
yason.py input.yml | faunus
~~~

produces an output file, `out.json`, with move statistics, system properties etc.
The script `yason.py` merely converts from YAML to JSON as the former, easier to read,
is used in all examples.

---

~~~ yaml
temperature: 300
geometry: { length: 100 }
mcloop: { macro: 10, micro: 1000 }
energy:
    - nonbonded_coulomblj:
        lennardjones: { mixing: LB }
        coulomb: { type: plain, epsr: 80, cutoff: 50 }
atomlist:
    - Na: { q:  1, r: 2, eps: 0.1, dp: 20 }
    - Cl: { q: -1, r: 2, eps: 0.1, dp: 20 }
moleculelist:
    - salt: { Ninit: 50, atoms: [Na,Cl], atomic: yes }
moves:
    - transrot: { molecule: salt }
~~~
