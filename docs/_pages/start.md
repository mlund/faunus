---
permalink: /start/
sidebar:
    nav: "docs"
---
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

# Getting Started

## Installation

### Requirements

- CMake 3.6 or higher
- C++14 compiler (clang, gcc etc.)
- Python 3.6 or higher with `ruamel_yaml` or `yaml`

### Compiling

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

<a name="input-output"></a>
## Input/Output 

Natively, Faunus input and output are [JSON formatted](http://json.org/example.html):

~~~ json
{ "atomlist": [
        { "Na+": { "q": 1.0, "mw": 22.99 } }
    ]
}
~~~

However, via the small helper script `yason.py`, JSON can be converted to/from
[YAML](http://www.yaml.org) which is less verbose, more readable and therefore
used throughout the documentation:

~~~ yaml
atomlist:
    - Na+: { q: 1.0, mw: 22.99 }
~~~

### Post-Processing 

Output (JSON) can be conveniently converted to
syntax highlighted YAML for better readability:

~~~ bash
yason.py --color out.json
~~~

For further processing of output or input, JSON (and YAML) can be read by
most programming languages. For example in python:

~~~ python
import json
with open('out.json') as f:
    d = json.load(f) # --> dict
    print( d['atomlist'][0]["Na+"]["mw"] ) # --> 22.99
~~~

<a name="running"></a>
## Running Faunus

Input is read from `stdin` and can either be JSON or,
via `yason.py`, also YAML:

~~~ bash
faunus < input.json # from json
yason.py minimal.yml | faunus # from yaml
~~~

