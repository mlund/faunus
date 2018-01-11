---
permalink: /input/
sidebar:
    nav: "docs"
---
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

# Handling Input and Output

In this section we describe the type of input you can use with Faunus.

## JSON Files

Faunus natively uses the [JSON format](http://json.org/example.html) for input and,
in many cases, also for output.
It is a structured format that can be easily loaded/saved from most programming
languages such as Python, C++ etc. Here's an example of an input file
for melted NaCl for the `bulk.cpp` program:

~~~ json
{
  "atomlist" : {
    "Na": {"q": 1, "sigma":3.3, "eps":0.3, "dp":1.0},
    "Cl": {"q":-1, "sigma":4.4, "eps":0.4, "dp":1.0}
  },
  "moleculelist" : { 
    "salt" : { "atoms":"Na Cl", "atomic":true, "Ninit":512}
  },
  "moves" : {
    "atomtranslate" : { "salt" : { "peratom":true } }
  },
  "energy" : {
    "nonbonded" : { "epsr": 80.0, "cutoff":14 }
  },
  "analysis" : {
    "xtcfile" : { "file": "traj.xtc", "nstep":20 },
    "pqrfile" : { "file": "confout.pqr" },
    "atomrdf" : { "nstep":10, "pairs" :
      [
        {"name1":"Na", "name2":"Cl", "dr":0.1, "file":"rdf1.dat"},
        {"name1":"Na", "name2":"Na", "dr":0.1, "file":"rdf2.dat"}
      ]
    }
  },
  "system" : {
      "temperature" : 1100,
      "geometry" : { "length":40 }
  }
}
~~~

### Reading From Python

~~~ python
import json
with open('move_out.json') as f:
    data = json.load(f)
    print( data['moves']['Single Particle Translation']['trials'] ) # --> 602580
~~~

### Reading From C++11

~~~ cpp
#include <iostream>
#include <json.hpp> // https://github.com/nlohmann/json

int main() {
    nlohmann::json j;
    j << std::cin;
    std::cout << j["moves"]["Single Particle Translation"]["trials"]; // --> 602580
}
~~~

### Syntax Highlighting

Python dictionaries can be directly converted to/from JSON which
has the advantage of providing syntax checking. If you prefer to
edit JSON files directly, simple syntax highlighting in the VIM editor
is possible by by adding 
the following to `~/.vimrc`:

    au! BufRead,BufNewFile *.json set filetype=javascript

### Conversion to/from YAML

JSON is a subset of the more general [YAML format](http://www.yaml.org) which is more human readablee.
To convert between JSON and YAML we
provide a script, `scripts/yason.py`, that can pass
YAML (and JSON) input files,

~~~ bash
yason.py < input.yml | ./myprogram
~~~

where we assume that `myprogram` reads JSON from standard input, which
in C++ would be achieved like this:

~~~ cpp
nlohmann::json j;
std::cin >> j;
~~~

The script can also be used oppositely to pretty print JSON as YAML,

~~~ bash
yason.py --color move_out.json
~~~

Example YAML file from a Monte Carlo move class:

~~~ yaml
moves:
  Single Particle Translation:
    acceptance: 0.63720
    atoms:
      Na: {acceptance: 63.7, dp: 1, mean displacement: 0.385}
    dir: [1, 1, 1]
    moves/particle: 20086
    relative time: 0.483
    trials: 602580
~~~

which is less verbose than the corresponding JSON file,

~~~ json
{
  "moves" {
    "Single Particle Translation": {
      "acceptance": 0.637,
      "atoms": {
        "Na": {
          "acceptance": 63.72,
          "dp": 1,
          "mean displacement": 0.385}
      },
      "dir": [ 1, 1, 1 ],
      "moves/particle": 20086,
      "relative time": 0.483,
      "trials": 602580
    },
  }
}
~~~


# Misc

<iframe src="https://docs.google.com/file/d/0BzpLUBrTxmurRzN0RnFZc2lhZFE/preview" width="640" height="385"></iframe>

## Make a new Faunus program

We'll here show three ways to link a new, external program to the Faunus library.

## Using CMake

Start by making a new directory (anywhere you want), put your source file there,
and make a `CMakeLists.txt` file telling CMake about the new executable.
For example:

~~~ bash
$ cd $HOME/newproject
$ cat hello.cpp

#include <faunus/faunus.h>
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

Note that all options selected when configuring faunus will be applied to `hello` as well,
and any changes to the faunus code base will trigger re-compilation upon running `make`.

# Coding Guidelines

## Naming Style

Object                          | Example                     | Comment
:------------------------------ | :-------------------------- | :----------
Types and namespaces            | `class AtomicTranslation;`  | Mixed case starting with upper case
Member functions, functions     | `double getTemperature();`  | Mixed case starting with lower case
Public variables (*avoid!*)     | `int numberOfParticles;`    | Mixed case starting with lower case
Private variables               | `int _i;`                   | Underscore prefix 

## Editing code

- Document code using Doxygen tags - the
  [Markdown](http://www.stack.nl/~dimitri/doxygen/markdown.html) syntax is recommended.
- Set your editor to use a *white space* indentation of *two*.
  - VIM: add the following to `.vimrc`:
  ~~~
  set expandtab
  set shiftwidth=4
  set softtabstop=4
  ~~~

## Design

Some books on C++ design,

- [C++ Coding Standards](http://en.wikipedia.org/wiki/Special:BookSources/0321113586)
  by Sutter and Alexandrescu as well as
- [Effective C++](http://en.wikipedia.org/wiki/Special:BookSources/0321334876)
  by Meyers.

A few basic tips,

- Generously use C++'s [assert()](http://www.cplusplus.com/reference/clibrary/cassert/assert)
  command to ease debugging
- Recycle code with polymorphic designs
- Use compile-time polymophism (templates) for speed limiting steps
- Use STL and Eigen where possible
- Hide data and functions as much as possible (i.e. make them private)
- Stride for [const-correctness](http://en.wikipedia.org/wiki/Const-correctness)
- Treat compiler warnings as errors

## Committing code

To contribute to the Faunus project, first make a *fork* of the
repository and request changes to be added via a *pull request*.
Further instructions can be found on GitHub:

- [Forking a repository](http://help.github.com/articles/fork-a-repo)
- [Pull requests](http://help.github.com/articles/using-pull-requests)

Before submitting changes, please make sure nothing is broken:

~~~ bash
make all
make manual
make test
~~~
