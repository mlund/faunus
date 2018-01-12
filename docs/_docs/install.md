---
---
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

# Installation

Faunus has been tested developed on MacOS and Linux but should work
on most unix operating systems.

## Requirements

- CMake 3.9 or higher
- C/C++14 compiler (clang, gcc etc.)

Optional packaged:

- Python 3.6 or higher with `ruamel_yaml` or `yaml`
- Message Parsing Interface (MPI)

## Compiling

Download the [latest release](https://github.com/mlund/faunus/releases/latest)
and perform the following steps in a terminal.
A set of dependencies will automatically be downloaded.

~~~ bash
cd faunus
cmake . [OPTIONS]
make faunus
~~~

The following options are available:

Option                             | Description
---------------------------------  | ---------------------------------------
`-DENABLE_MPI=OFF`                 | Enable MPI
`-DENABLE_OPENMP=OFF`              | Enable OpenMP support
`-DENABLE_PYTHON=ON`               | Build python bindings (experimental)
`-DENABLE_POWERSASA=OFF`           | Enable SASA routines (external download)
`-DCMAKE_BUILD_TYPE=RelWithDebInfo`| Alternatives: `Debug` or `Release` (faster)
`-DCMAKE_CXX_FLAGS_RELEASE="..."`  | Compiler options for Release mode
`-DCMAKE_CXX_FLAGS_DEBUG="..."`    | Compiler options for Debug mode
`-DMYPLAYGROUND="absolute path"`   | Add additional source directory

## Libraries in odd locations

Should you have multiple compilers or python distributions, be specific:

~~~ bash
CC=clang CXX=clang++ cmake . \
  -DPYTHON_INCLUDE_DIR=/opt/include/python3.6 \
  -DPYTHON_LIBRARY=/opt/lib/libpython3.6.dylib
~~~
