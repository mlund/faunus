---
---
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

# Installation

## Requirements

- CMake 3.9 or higher
- C++14 compiler (clang, gcc etc.)

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

Should you have multiple compilers or python distributions, be specific:

~~~ bash
CXX=clang++ cmake . -DPYTHON_INCLUDE_DIR=/opt/include/python3.6 -DPYTHON_LIBRARY=/opt/lib/libpython3.6.dylib
~~~


