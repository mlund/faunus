[![Documentation](https://img.shields.io/badge/documentation-red.svg)](http://mlund.github.io/faunus/)
[![Anaconda-Server Badge](https://anaconda.org/teokem/faunus/badges/installer/conda.svg)](https://conda.anaconda.org/teokem/)
[![Build Status](https://travis-ci.org/mlund/faunus.svg?branch=master)](https://travis-ci.org/mlund/faunus)
[![License: MIT](https://img.shields.io/badge/License-MIT-red.svg)](https://opensource.org/licenses/MIT)

Welcome to Faunus
=================

Faunus is a C++ framework for Metropolis Monte Carlo simulations of
molecular systems. Below is a brief overview of features:

- Canonical, Grand Canonical, Isobaric-Isothermal statistical mechanical ensembles
- General hamiltonian **parallel tempering** (temperature, screening length, bonds etc.)
- **Anisotropic** atoms (multipoles, sphero-cylinders, capped particles)
- Speciation moves (pKa prediction, Hofmeister effects etc.)
- Highy modular
- Free and open source

Installing
===========

On mac or linux, install using [conda](https://conda.io/miniconda.html):

    $ conda install -c teokem faunus

Documentation
=============

http://mlund.github.io/faunus

Licence
=======

Copyright 2002-2018 Mikael Lund

Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software
is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
