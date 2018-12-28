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
to use, and to modify. The code is written in C++ and Python bindings
are available.
The development is a team effort with, in reverse chronological order,
many valiant contributions from:

_Vidar Aspelin, Stefan Hervø-Hansen,
Richard Chudoba, Niels Kouwenhoven,
Coralie Pasquier, Lukáš Sukeník,
Giulio Tesei, Alexei Abrikossov,
João Henriques, Björn Stenqvist,
Axel Thuresson, Robert Vácha,
Magnus Ullner, Chris Evers,
Anıl Kurut, André Teixeira,
Christophe Labbez, Ondrej Marsalek,
Martin Trulsson, Björn Persson,
[Mikael Lund](http://www.teokem.lu.se/~mikael)_

Should you find Faunus useful, please consider supporting us by crediting:

- Stenqvist _et al._ [_Molecular Simulation 2013, 39:1233_](http://dx.doi.org/10/nvn)
- Lund, M. _et al._ [_Source Code Biol. Med., 2008, 3:1_](http://dx.doi.org/10/dfqgch)

## Quick Start

Simulations are set up using YAML or JSON files, see for example
[minimal.yml](https://github.com/mlund/faunus/blob/master/examples/minimal.yml),
for a Metropolis Monte Carlo simulation of charged Lennard-Jones particles in a cubic PBC box.
Running with

~~~ bash
yason.py minimal.yml | faunus
~~~

produces an output file, `out.json`, with move statistics, system properties etc.
The script `yason.py` merely converts from YAML to JSON as the former, easier to read,
is used in all examples.
For more examples, see the `examples` folder.

To open the user-guide in a browser, type:

~~~ bash
faunus-manual
~~~

