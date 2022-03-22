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

Starting from version 2.1, we adhere to [semantic versioning](https://semver.org).
**Note:** Updating to a newer version is often delayed, and if the version you're after
is not on Conda, consider an alternative method below.


## Using Docker

We provide a [`Dockerfile`](https://github.com/mlund/faunus/blob/master/scripts/Dockerfile)
that builds the main branch in a [Jupyter](https://jupyter.org) environment:

~~~ bash
curl -s https://raw.githubusercontent.com/mlund/faunus/master/scripts/Dockerfile | docker build -t faunus -
docker run -it -p 8888:8888 faunus # paste generated url into webbrowser
~~~


## Building from source code

Faunus is continuously [tested](https://travis-ci.org/mlund/faunus) on macOS/Linux,
but should compile on most unix operating systems and possibly under Cygwin (Windows).

### Requirements

- CMake 3.16+
- C/C++17 compiler (Clang 5+, GCC 7+, etc.)
- Python 3.6+ with the following packages:
  - `jinja2`, `ruamel_yaml` or `yaml`

The following are optional:

- `jsonschema` (for validating input - highly recommended)
- `pandoc` (for building manual)
- `pypandoc` (for building manual)
- `BeautifulSoup4` (for building manual)
- Message Passing Interface (MPI)

**macOS tip:**
Apple's developer tools, Xcode, include clang;
CMake can be installed with an
[Installer package](https://cmake.org/download) from Kitware, or using
[Homebrew](https://brew.sh), or
[(mini)conda](https://conda.io/docs/user-guide/install/index.html)

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
`-DENABLE_TESTS=ON`                  | Enable unittesting
`-DENABLE_PYTHON=OFF`                | Build python bindings (experimental)
`-DENABLE_FREESASA=ON`               | Enable SASA routines (external download)
`-DENABLE_TBB=OFF`                   | Build with Intel Threading Building Blocks (experimental)
`-DENABLE_PCG=OFF`                   | Use PCG random number generator instead of C++'s Mersenne Twister
`-DBUILD_STATIC=OFF`                 | Build statically linked binaries
`-DCMAKE_BUILD_TYPE=Release`         | Alternatives: `Debug` or `RelWithDebInfo`
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

In addition to pandoc, a TeX Live installation containing XeLaTeX is required to build the PDF manual.
The manual is supposed to be typeset with
[EB Garamond](https://github.com/octaviopardo/EBGaramond12/tree/master/fonts/otf),
[Garamond Math](https://github.com/YuanshengZhao/Garamond-Math/blob/master/Garamond-Math.otf) and
[Fira Code](https://github.com/tonsky/FiraCode/releases/download/2/FiraCode_2.zip)
fonts thus they have to be available in your system. Alternatively, you can tweak the font options
in the `header.md` file.

~~~ bash
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
rm -fR CMakeCache.txt CMakeFiles _deps
~~~

### Intel Threading Building Blocks

If `ENABLE_TBB=on`, TBB may be used for threaded simulations which may or may not be
advantageous, depending on the system.
By default, an unspecified and possibly outdated version of TBB will be downloaded and build.
Alternatively you can use an existing installation _via_ `TBB_DIR`:

~~~ bash
cmake -DENABLE_TBB=on -DTBB_DIR={tbb-root}/lib/cmake/TBB
~~~

where `{tbb-root}` is the installation directory of TBB, _e.g._ `/usr/local`.

# Development

The development of Faunus is done mainly in Jetbrain's [CLion](https://www.jetbrains.com/clion)
(free academic license) but any other IDE or merely a text editor can be used.
We recommend to use tools that respect the provided `.clang-format` which will ease merging
changes into the codebase, see below.

## Code Style

If you plan to contribute to Faunus it is recommended to activate the
pre-commit hook for automatic styling of all changes:

``` bash
cd faunus
./scripts/git-pre-commit-format install
```

This requires `clang-format` which may also be directly used in IDEs
such as CLion. In the top-level directory of Faunus you will find
the style configuration file [`.clang-format`](https://github.com/mlund/faunus/blob/master/.clang-format)

Also, adhere to the following naming conventions:

Style        | Elements
------------ | -------------------------
`CamelCase`  | classes, namespaces
`camelBack`  | functions
`snake_case` | variables


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
