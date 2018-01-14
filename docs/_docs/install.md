---
---
# Installation

Faunus has been tested developed on MacOS/Linux, but should work
on any unix operating systems and _possibly_ under Cygwin (Windows).

## Requirements

- CMake 3.9 or higher
- C/C++14 compiler (clang, gcc etc.)

Optional packages:

- Python 3.6 or higher with `ruamel_yaml` or `yaml`
- Message Passing Interface (MPI)

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

## Linking an external program

Start by making a new directory (anywhere you want), put your .cpp source file there,
and make a `CMakeLists.txt` file telling CMake about the new executable.
For example:

~~~ bash
$ cd $HOME/newproject
$ cat hello.cpp

#include "faunus.h"
int main() {
  Faunus::Point a(0,0,0);
}

$ echo 'fau_example(hello "./" hello.cpp)' > CMakeLists.txt
~~~

Return to the main faunus directory and rerun `cmake` with the following command:

~~~ bash
cd $HOME/faunus
cmake . -DMYPLAYGROUND=$HOME/newproject  # absolute path!
~~~

That's it! A `Makefile` for your new target, `hello`, has been generated and you can compile
directly from the `newproject` directory:

~~~ bash
cd $HOME/newproject
make
~~~

All options selected when configuring faunus will be applied to `hello` as well,
and changes to the faunus code base will trigger re-compilation upon running `make`.

