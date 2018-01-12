---
---
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

# Installation

Faunus has been tested developed on MacOS/Linux, but should work
on any unix operating systems and _possibly_ under Cygwin (Windows).

## Requirements

- CMake 3.9 or higher
- C/C++14 compiler (clang, gcc etc.)

Optional packages:

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
`-DENABLE_OPENMP=ON`               | Enable OpenMP support
`-DENABLE_PYTHON=ON`               | Build python bindings (experimental)
`-DENABLE_POWERSASA=ON`            | Enable SASA routines (external download)
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

Another example for compiling with Intel C++ in _Release_ mode (faster, less assertions):

~~~ bash
CXX=icpc CC=icc cmake . -DCMAKE_BUILD_TYPE=Release
make
~~~

If you experience python issues on macOS, the linked python library can be probed and,
if needed, renamed:

~~~ bash
otool -L pyfaunus.so
install_name_tool -change libpython3.6.dylib $HOME/miniconda/lib/libpython3.6.dylib pyfaunus.so
~~~

## Resetting the build system

If you need to change the compiler or for another reason want to reset the build system, do:

~~~ bash
make clean
rm -fR CMakeCache.txt CMakeFiles
~~~

