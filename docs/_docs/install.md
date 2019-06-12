# Installing

## Using Conda

For macOS and Linux x86-64, precompiled binary packages are available
via [(mini)conda](https://conda.io/docs/user-guide/install/index.html):

~~~ bash
conda config --add channels conda-forge
conda install faunus
~~~

In addition to the `faunus` executable, this installs a set of examples in `share/faunus`,
as well as python bindings.
To _update_ an existing installation, use

~~~ bash
faunus --version      # show version string
conda search faunus   # show (new) revisions
conda upgrade faunus
~~~

For the adventurous, the latest development version is available with `conda install -c teokem faunus`. 

Starting from version 2.1, we adhere to [semantic versioning](https://semver.org).

## Building from source code

Faunus is continuously [tested](https://travis-ci.org/mlund/faunus) on macOS/Linux,
but should compile on most unix operating systems and possibly under Cygwin (Windows).

### Requirements

- CMake 3.11+
- C/C++14 compiler (Clang 3.5+, GCC 6+, etc.)
- Python 3.6+ with the following packages:
  - `ruamel_yaml` or `yaml`

The following are optional:

- `jinja2`
- `pandoc`
- `pypandoc`
- `BeautifulSoup4`
- Message Passing Interface (MPI)

**macOS tip:**
Apple's developer tools, Xcode, include clang and
CMake can be installed with an
[Installer package](https://cmake.org/download) from Kitware, or using
[MacPorts](http://www.macports.org),
[Homebrew](https://brew.sh), or
[(mini)conda](https://conda.io/docs/user-guide/install/index.html)
{: .notice--info}

### Compiling

Download the [latest release](https://github.com/mlund/faunus/releases/latest)
or [the developer branch](https://github.com/mlund/faunus/archive/master.zip)
and build using cmake:

~~~ bash
cd faunus
cmake . [OPTIONS]
make faunus
make usagetips # requires `pandoc`, `pypandoc`, `BeautifulSoup4`
~~~

Use `make help` to see all build targets.

The following options are available:

CMake Option                         | Description
------------------------------------ | ---------------------------------------
`-DENABLE_MPI=OFF`                   | Enable MPI
`-DENABLE_OPENMP=ON`                 | Enable OpenMP support
`-DENABLE_PYTHON=ON`                 | Build python bindings (experimental)
`-DENABLE_POWERSASA=ON`              | Enable SASA routines (external download)
`-DBUILD_STATIC=OFF`                 | Build statically linked binaries
`-DCMAKE_BUILD_TYPE=RelWithDebInfo`  | Alternatives: `Debug` or `Release` (faster)
`-DCMAKE_CXX_FLAGS_RELEASE="..."`    | Compiler options for Release mode
`-DCMAKE_CXX_FLAGS_DEBUG="..."`      | Compiler options for Debug mode
`-DCMAKE_INSTALL_PREFIX:PATH="..."`  | Install location (default: `/usr/local`)
`-DPYTHON_EXECUTABLE="..."`          | Full path to Python executable
`-DPYTHON_INCLUDE_DIR="..."`         | Full path to python headers
`-DPYTHON_LIBRARY="..."`             | Full path to python library, i.e. libpythonX.dylib/so

### Compiling the Manual

Pandoc is required to build the HTML manual: 

~~~ bash
make manual_html
~~~

In addition to pandoc, a TeX Live installation is required to build
the PDF manual. Garamond fonts must be available:

~~~ bash
wget https://tug.org/fonts/getnonfreefonts/install-getnonfreefonts
sudo texlua install-getnonfreefonts
sudo getnonfreefonts garamond --sys
make manual
~~~

### Python libraries in odd locations

Should there be multiple compilers or python distributions, be specific:

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

# Development

## Code Style

If you plan to contribute to Faunus it is recommended to activate the
pre-commit hook for automatic styling of all changes:

``` bash
cd faunus
./scripts/git-pre-commit-format install
```

This requires `clang-format` which may also be directly used in IDE's
such as CLion. In the top-level directory of Faunus you will find
the style configuration file `.clang-format`


## Creating a conda package (development usage)

The basic steps for creating a conda package is outlined below, albeit
details depend on the build environment. See also the `.travis.yml`
configuration file in the main repository.

~~~ bash
conda config --add channels conda-forge
conda install conda-build anaconda-client
cd scripts/
conda-build .
anaconda login
anaconda upload -u USER ... # see output from build step
~~~

Instead of uploading to anaconda.org, install a local copy directly after the build step above:

~~~ bash
conda install -c USER faunus --use-local
~~~

