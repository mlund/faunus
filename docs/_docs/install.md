---
---
[![Edit](https://img.shields.io/badge/Github-Improve_this_page-orange.svg)]({{site.github.repository_url}}/blob/master/docs/{{page.path}})

# Installing

## Using Conda

For macOS and Linux x86_64, precompiled binary packages are available
via [(mini)conda](https://conda.io/docs/user-guide/install/index.html):

~~~ bash
conda config --add channels conda-forge
conda install -c teokem faunus
~~~

In addition to the `faunus` program, this installs a set of examples in `share/faunus`.

## Building from source code

Faunus is continuously [tested](https://travis-ci.org/mlund/neofaunus) on macOS/Linux,
but should compile on most unix operating systems and possibly under Cygwin (Windows).

### Requirements

- CMake 3.9+
- C/C++14 compiler (Clang 3.9+, GCC 6+, etc.)

The following are optional:

- Python 3.6+ with `ruamel_yaml` or `yaml`
- Message Passing Interface (MPI)
- Pandoc (for building documentation)

**macOS tip:**
Apple's developer tools, Xcode, is an easy way obtain
clang on macOS. CMake can be installed with an
[Installer package](https://cmake.org/download) from Kitware, or using
[MacPorts](http://www.macports.org), or
[Homebrew](https://brew.sh).
{: .notice--info}

### Compiling

Download the [latest release](https://github.com/mlund/faunus/releases/latest)
and perform the following steps in a terminal.
A set of dependencies will automatically be downloaded.

~~~ bash
cd faunus
cmake . [OPTIONS]
make
make tests
make install
~~~

The following options are available:

Option                               | Description
------------------------------------ | ---------------------------------------
`-DENABLE_MPI=OFF`                   | Enable MPI
`-DENABLE_OPENMP=OFF`                | Enable OpenMP support
`-DENABLE_PYTHON=ON`                 | Build python bindings (experimental)
`-DENABLE_POWERSASA=ON`              | Enable SASA routines (external download)
`-DCMAKE_BUILD_TYPE=RelWithDebInfo`  | Alternatives: `Debug` or `Release` (faster)
`-DCMAKE_CXX_FLAGS_RELEASE="..."`    | Compiler options for Release mode
`-DCMAKE_CXX_FLAGS_DEBUG="..."`      | Compiler options for Debug mode
`-DCMAKE_INSTALL_PREFIX:PATH="..."`  | Install location (default: /usr/local)
`-DMYPLAYGROUND="absolute path"`     | Add additional source directory

### Libraries in odd locations

Should you have multiple compilers or python distributions, be specific:

~~~ bash
CC=clang CXX=clang++ cmake . \
  -DPYTHON_INCLUDE_DIR=/opt/include/python3.6 \
  -DPYTHON_LIBRARY=/opt/lib/libpython3.6.dylib
~~~

Another example for compiling with Intel C++ in _Release_ mode (faster, less assertions):

~~~ bash
CXX=icpc CC=icc cmake . -DCMAKE_BUILD_TYPE=Release
~~~

For solving python issues on macOS, the linked python library can be probed and,
if needed, renamed:

~~~ bash
otool -L pyfaunus.so
install_name_tool -change libpython3.6.dylib \
  $HOME/miniconda/lib/libpython3.6.dylib pyfaunus.so
~~~

### Resetting the build system

If you need to change the compiler or for another reason want to reset the build system, do:

~~~ bash
make clean
rm -fR CMakeCache.txt CMakeFiles
~~~

## Creating a conda package (expert usage)

To create a precompiled package for Anaconda, do the following steps:

~~~ bash
conda config --add channels conda-forge
conda install conda-build anaconda-client
cd scripts/
conda-build .
anaconda login
anaconda upload -u USER ... # see output from build step
~~~

Alternatively, instead of uploading to anaconda.org you
can install a local copy directly after the build step above:

~~~ bash
conda install -c USER faunus --use-local
~~~

### Requirements

Building a conda package require a full LaTeX installation in order to create
the PDF manual. Do ensure that the Garamond fonts are available:

~~~ bash
wget https://tug.org/fonts/getnonfreefonts/install-getnonfreefonts
sudo texlua install-getnonfreefonts
sudo getnonfreefonts garamond --sys
~~~
