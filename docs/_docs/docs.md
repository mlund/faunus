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
to use and to modify. The code is written in C++ and (limited) Python bindings
are available.
The development is a **team effort** with, in chronological order,
many valiant contributions from:

_[Mikael Lund](http://www.teokem.lu.se/~mikael),
Björn Persson, Martin Trulsson,
Ondrej Marsalek, Christophe Labbez, André Teixeira,
Anıl Kurut, Chris Evers, Magnus Ullner,
Robert Vácha, Axel Thuresson, Björn Stenqvist,
João Henriques, Alexei Abrikossov, Giulio Tesei,
Lukáš Sukeník, Coralie Pasquier, Niels Kouwenhoven,
Richard Chudoba, Stefan Hervø-Hansen_.

### Supporting Faunus

Should you find Faunus useful, do consider supporting us by crediting:

- Stenqvist _et al._ [_Molecular Simulation 2013, 39:1233_](http://dx.doi.org/10/nvn)
- Lund, M. _et al._ [_Source Code Biol. Med., 2008, 3:1_](http://dx.doi.org/10/dfqgch)

## Quick Start

Simulations are set up using YAML or JSON files, see for example
[minimal.yml](https://github.com/mlund/faunus/blob/master/examples/minimal.yml),
for a Metropolis Monte Carlo simulation
of charged Lennard-Jones particles in a cubic PBC box. Running is with

~~~ bash
yason.py input.yml | faunus
~~~

produces an output file, `out.json`, with move statistics, system properties etc.
The script `yason.py` merely converts from YAML to JSON as the former, easier to read,
is used in all examples.
For more examples, check out the `examples/` folder.

