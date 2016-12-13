Welcome to Faunus
=================

Faunus is a C++ framework for Metropolis Monte Carlo simulations of
molecular systems. Below is a brief overview of features:

- Canonical, Grand Canonical, Isobaric-Isothermal statistical mechanical ensembles
- General hamiltonian **parallel tempering** (temperature, screening length, bonds etc.)
- **Anisotropic** atoms (multipoles, sphero-cylinders)
- Ion titration moves (pKa prediction, Hofmeister effects etc.)
- Highy modular
- Free and open source (**GPL**)

Detailed information and support:

- Manual: <http://faunus.sourceforge.net>
- Support: <http://github.com/mlund/faunus/issues>
- Source code: <http://github.com/mlund/faunus>

Get the latest revision from GitHub:

    $ git clone https://github.com/mlund/faunus.git

Requirements
============

- C/C++11 compiler (clang3.5+, gcc4.9+, intel14+, ...)
- CMake 2.8.5+

Optional:

- Doxygen (for code manual)
- MPI (for parallelisation)
- Python (for python module, experimental)

Developed and tested on Linux and MacOS X.

Compiling
=========

    $ cmake . [options]
    $ make [help]
    $ make test (or use 'ctest -V' for verbose test output)

Build options
-------------

Option                             | Description
:--------------------------------- | :----------------------------------------
`-DENABLE_MPI=OFF`                 | Build MPI programs (parallel tempering etc.)
`-DENABLE_OPENMP=OFF`              | Enable OpenMP support
`-DENABLE_STATIC=OFF`              | Static linkage of faunus as opposed to default dynamic linkage
`-DENABLE_UNICODE=ON`              | Use Unicode UTF-16 encoding for pretty output
`-DENABLE_PYTHON=ON`               | Build python bindings (experimental)
`-DENABLE_POWERSASA=OFF`           | Enable SASA routines (external download)
`-DCMAKE_BUILD_TYPE=RelWithDebInfo`| Alternatives: `Debug` or `Release` (faster)
`-DCMAKE_CXX_FLAGS_RELEASE="..."`  | Compiler options for Release mode
`-DCMAKE_CXX_FLAGS_DEBUG="..."`    | Compiler options for Debug mode
`-DMYPLAYGROUND="absolute path"`   | Add additional source directory

Example: Intel's C++ compiler
-----------------------------

    $ CXX=icpc CC=icc cmake . -DCMAKE_BUILD_TYPE=Release
    $ make

Example: Libraries in odd locations
-----------------------------------

    $ LDFLAGS=-L/sw/lib CPPFLAGS=-I/sw/include cmake .

Examples: Python support
------------------------
Faunus has preliminary python support -- see this Jupyter [notebook](scripts/pyfaunus-test.ipynb)
for examples.
If several python distributions are installed, the build system may be guided to specific
libraries using i.e.:

    $ cmake . -DPYTHON_INCLUDE_DIR=$HOME/miniconda/include/python2.7 -DPYTHON_LIBRARY=$HOME/minoconda/lib/libpython2.7.dylib
    $ make pyfaunus

On OSX the linked python library can be probed and, if needed, renamed like this:

    $ otool -L pyfaunus.so
    $ install_name_tool -change libpython2.7.dylib $HOME/miniconda/lib/libpython2.7.dylib pyfaunus.so

Resetting the build system
--------------------------

    $ make clean
    $ rm -fR CMakeCache.txt CMakeFiles

Contributors
============

In chronological order:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Mikael Lund         Bjorn Persson       Martin Trulsson    
Ondrej Marsalek     Christophe Labbez     Andre Teixeira     
  Anil Kurut           Chris Evers         Magnus Ullner      
 Robert Vacha         Axel Thuresson      Bjorn Stenqvist
 Joao Henriques     Alexei Abrikossov      Giulio Tesei
 Lukas Sukenik      Niels Kouwenhoven
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Citing Faunus
=============

- Stenqvist et al.
  _Molecular Simulation 2013, 39:1233_

  [![DOI](https://img.shields.io/badge/DOI-10%2Fnvn-orange.svg)](http://dx.doi.org/10/nvn)

- Lund, M., Persson, B., Trulsson, M.
  _Source Code Biol. Med., 2008, 3:1_

  [![DOI](https://img.shields.io/badge/DOI-10%2Fdfqgch-orange.svg)](http://dx.doi.org/10/dfqgch)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 Faunus - A Framework for Molecular Modelling 
 Copyright (C) 2002-2016 Mikael Lund

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or 
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

