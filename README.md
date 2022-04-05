[![Documentation](https://readthedocs.org/projects/faunus/badge/?version=latest)](https://faunus.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://app.travis-ci.com/mlund/faunus.svg?branch=master)](https://app.travis-ci.com/mlund/faunus)
[![License: MIT](https://img.shields.io/badge/License-MIT-brightgreen.svg)](https://opensource.org/licenses/MIT)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/faunus/badges/installer/conda.svg)](https://anaconda.org/conda-forge/faunus)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/2ac7217d373a4f34a2dae2d912c9d1a1)](https://www.codacy.com/app/mlund/faunus?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=mlund/faunus&amp;utm_campaign=Badge_Grade)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5235137.svg)](https://doi.org/10.5281/zenodo.5235137)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/faunus/badges/platforms.svg)](https://anaconda.org/conda-forge/faunus)


Welcome to Faunus
=================

Faunus is a molecular simulation package for Metropolis Monte Carlo simulations of
molecular systems. Below is a brief overview of features:

- Canonical, Grand Canonical, Isobaric-Isothermal statistical mechanical ensembles
- General hamiltonian parallel tempering (temperature, screening length, bonds etc.)
- Anisotropic atoms (multipoles, sphero-cylinders, capped particles)
- Speciation moves (pKa prediction, Hofmeister effects, general equilibrium reactions etc.)
- Parallelise flat histogram method for sampling free energies along 1D/2D reaction coordinates
- Highy modular
- Free and open source

Installing
===========

On macOS or Linux, install latest release using [conda](https://conda.io/miniconda.html):

    conda install -c conda-forge faunus

Or [build development branch from source](https://faunus.readthedocs.io/en/latest/_docs/install.html#building-from-source-code).

Documentation
=============

- [latest release](https://github.com/mlund/faunus/releases/latest)
- [development branch](https://faunus.readthedocs.io/en/latest/?badge=latest)
- installed version: type `faunus-manual`

Licence
=======

Copyright 2002-2022 Mikael Lund

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
