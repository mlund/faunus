Welcome to Faunus
=================

Faunus is a C++ framework for Metropolis Monte Carlo simulations of
molecular systems. Below is a brief overview of features:

- Canonical, Grand Canonical, Isobaric-Isothermal statistical mechanical ensembles
- General hamiltonian **parallel tempering** (temperature, screening length, bonds etc.)
- 4D simulations on **hyperspheres**
- **Anisotropic** atoms (multipoles, sphero-cylinders)
- Ion titration moves (pKa prediction, Hofmeister effects etc.)
- Highly expandable!
- Free and open source (**GPL**)

Detailed information and support:

- Website: <http://github.com/mlund/faunus>
- Documentation: <http://faunus.sourceforge.net/doxyhtml>
- Support: <http://github.com/mlund/faunus/issues>

Get the latest resivions from GitHub:

    $ git clone https://github.com/mlund/faunus.git

Requirements
============

- C/C++11 compiler (clang3+, gcc4.6+, intel13+, ...)
- CMake 2.8+

Optional:

- OpenBabel 2 (numerous molecular formats)
- Doxygen and Graphviz (for code manual)
- Xcode (for neat environment on macos x)

Developed and tested on Linux and MacOS X.

Compiling                                         {#compiling}
=========

    $ cmake . [options]
    $ make [help]
    $ make test (or use 'ctest -V' for verbose test output)

Build options
-------------

Option                             | Description
:--------------------------------- | :----------------------------------------
`-DENABLE_BABEL=OFF`               | Use OpenBabel for file I/O
`-DENABLE_MPI=OFF`                 | Build MPI programs (parallel tempering etc.)
`-DENABLE_OPENMP=OFF`              | Enable OpenMP support
`-DENABLE_STATIC=OFF`              | Static linkage of faunus as opposed to default dynamic linkage
`-DENABLE_TWISTER=OFF`             | Use Mersenne Twister for random numbers
`-DENABLE_UNICODE=ON`              | Use Unicode UTF-16 encoding for pretty output
`-DCMAKE_BUILD_TYPE=RelWithDebInfo`| Alternatives: `Debug` or `Release` (faster)
`-DCMAKE_CXX_FLAGS_RELEASE="..."`  | Compiler options for Release mode
`-DCMAKE_CXX_FLAGS_DEBUG="..."`    | Compiler options for Debug mode
`-DMYPLAYGROUND="..."`             | Add additional source directory (use absolute path)

Example: Intel's C++ compiler
-----------------------------

    $ CXX=icpc CC=icc cmake . -DCMAKE_BUILD_TYPE=Release
    $ make

Example: Libraries in odd locations
-----------------------------------

    $ LDFLAGS=-L/sw/lib CPPFLAGS=-I/sw/include/openbabel-2.0 cmake . -DENABLE_BABEL=on

Resetting the build system
--------------------------

    $ make clean
    $ rm CMakeCache.txt
    $ (rm -fR CMakeFiles/)

Contributors
============

In chronological order:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Mikael Lund         Bjorn Persson       Martin Trulsson    
Ondrej Marsalek     Christophe Labbez     Andre Teixeira     
  Anil Kurut           Chris Evers         Magnus Ullner      
 Robert Vacha         Axel Thuresson      Bjorn Stenqvist
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Citing Faunus
=============

- Lund, M., Persson, B., Trulsson, M.
  _Source Code Biol. Med., 2008, 3:1_
  [[DOI]](http://dx.doi.org/10/dfqgch)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 Faunus - A Framework for Molecular Modelling 
 Copyright (C) 2002-2013 Mikael Lund 

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


