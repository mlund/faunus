---
---
[![Edit](https://img.shields.io/badge/Github-Improve_this_page-orange.svg)]({{site.github.repository_url}}/blob/master/docs/{{page.path}})

# Installing

## Using Conda

For macOS and Linux x86-64, precompiled binary packages are available
via [(mini)conda](https://conda.io/docs/user-guide/install/index.html):

~~~ bash
conda config --add channels conda-forge
conda install -c teokem faunus
~~~

In addition to the `faunus` executable, this installs a set of examples in `share/faunus`,
as well as python bindings.
To _update_ an existing installation, use

~~~ bash
faunus --version               # show version string
conda search -c teokem faunus  # show (new) revisions
conda upgrade -c teokem faunus
~~~

Starting from version 2.1, we adhere to [semantic versioning](https://semver.org).

## Building from source code

Faunus is continuously [tested](https://travis-ci.org/mlund/faunus) on macOS/Linux,
but should compile on most unix operating systems and possibly under Cygwin (Windows).

### Requirements

- CMake 3.9+
- C/C++14 compiler (Clang 3.5+, GCC 6+, etc.)
- Python 3.6+ with `ruamel_yaml` or `yaml`

The following are optional:

- Message Passing Interface (MPI)
- Pandoc (for building documentation)

**macOS tip:**
Apple's developer tools, Xcode, is a quick way obtain
clang on macOS. CMake can be installed with an
[Installer package](https://cmake.org/download) from Kitware, or using
[MacPorts](http://www.macports.org),
[Homebrew](https://brew.sh), or
[(mini)conda](https://conda.io/docs/user-guide/install/index.html)
{: .notice--info}

### Compiling

Download the [latest release](https://github.com/mlund/faunus/releases/latest)
or [the developer branch](https://github.com/mlund/faunus/archive/master.zip)
and perform the following steps in a terminal.
Dependencies will automatically be downloaded.

~~~ bash
cd faunus
cmake . [OPTIONS]
make
make test
make install
~~~

The following options are available:

CMake Option                         | Description
------------------------------------ | ---------------------------------------
`-DENABLE_MPI=OFF`                   | Enable MPI
`-DENABLE_OPENMP=OFF`                | Enable OpenMP support
`-DENABLE_PYTHON=ON`                 | Build python bindings (experimental)
`-DENABLE_POWERSASA=ON`              | Enable SASA routines (external download)
`-DCMAKE_BUILD_TYPE=RelWithDebInfo`  | Alternatives: `Debug` or `Release` (faster)
`-DCMAKE_CXX_FLAGS_RELEASE="..."`    | Compiler options for Release mode
`-DCMAKE_CXX_FLAGS_DEBUG="..."`      | Compiler options for Debug mode
`-DCMAKE_INSTALL_PREFIX:PATH="..."`  | Install location (default: `/usr/local`)
`-DMYPLAYGROUND="absolute path"`     | Add additional source directory
`-DPYTHON_EXECUTABLE="..."`          | Full path to Python executable
`-DPYTHON_INCLUDE_DIR="..."`         | Full path to python headers
`-DPYTHON_LIBRARY="..."`             | Full path to python library, i.e. libpythonX.dylib/so


### Python libraries in odd locations

Should you have multiple compilers or python distributions, be specific:

~~~ bash
CC=/opt/bin/clang CXX=/opt/bin/clang++ cmake . \
  -DPYTHON_EXECUTABLE=/opt/bin/python3 \
  -DPYTHON_INCLUDE_DIR=/opt/include/python3.6 \
  -DPYTHON_LIBRARY=/opt/lib/libpython3.6.dylib
~~~

For solving python issues on macOS, the linked python library can be probed and,
if needed, renamed:

~~~ bash
otool -L pyfaunus.so
install_name_tool -change libpython3.6.dylib \
  $HOME/miniconda/lib/libpython3.6.dylib pyfaunus.so
~~~

### Resetting the build system

To change the compiler or for another reason reset the build system, do:

~~~ bash
make clean
rm -fR CMakeCache.txt CMakeFiles
~~~

## Creating a conda package (advanced usage)

The basic steps for creating a conda package is outlined below, albeit
details may depend on your build environment. See also the `.travis.yml`
configuration file in the main repository.

~~~ bash
conda config --add channels conda-forge
conda install conda-build anaconda-client
cd scripts/
conda-build .
anaconda login
anaconda upload -u USER ... # see output from build step
~~~

Instead of uploading to anaconda.org you
can install a local copy directly after the build step above:

~~~ bash
conda install -c USER faunus --use-local
~~~

### Requirements

Building a conda package may require a TeX Live installation to create
the PDF manual. Ensure that Garamond fonts are available:

~~~ bash
wget https://tug.org/fonts/getnonfreefonts/install-getnonfreefonts
sudo texlua install-getnonfreefonts
sudo getnonfreefonts garamond --sys
~~~
