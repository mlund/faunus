---
---
[![Edit](https://img.shields.io/badge/Github-Improve_this_page-orange.svg)]({{site.github.repository_url}}/blob/master/docs/{{page.path}})
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

# Running Simulations

The main program for running simulations is `faunus` and should be available
from the command line after installation. For general usage, type:

~~~ bash
faunus --help
~~~

Input is read either from `stdin` or from a JSON formatted file. Some examples:

~~~ bash
faunus < input.json              # input from stdin
faunus -i in.json -o out.json -q # file input/output and be quiet
~~~

Via the script `yason.py`, see below, [YAML](http://www.yaml.org)
formatted input can be passed:

~~~ bash
yason.py in.yml | faunus # from yaml
~~~

## Input and Output

Natively, input and output are [JSON formatted](http://json.org/example.html):

~~~ json
{ "atomlist": [
        { "Na+": { "q": 1.0, "mw": 22.99 } }
    ]
}
~~~

However, via the helper script `yason.py`, JSON can be converted to/from
[YAML](http://www.yaml.org) which is less verbose, more readable and therefore
used throughout the documentation:

~~~ yaml
atomlist:
    - Na+: { q: 1.0, mw: 22.99 }
~~~

### Post-Processing

JSON formatted output can conveniently be converted to
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

## Restarting

~~~ bash
faunus --input in.json --state state
~~~

## Message Passing Interface

Running with MPI spawns `nproc` processes that may or may not communicate
with each other. If `nproc>1`, input and output files are prefixed with
`mpi{rank}.` where `{rank}` is the rank or process number.

The following starts two processes, reading input from `mpi0.in.json` and
`mpi1.in.json`. All output files, including those from any analysis are
prefixed with `mpi0.` and `mpi1.`.

~~~
mpirun -np 2 ./faunus -i in.json
~~~

If all processes take the same input:

~~~
mpirun -np 2 ./faunus --nopfx --input in.json
mpirun -np 2 --stdin all ./faunus < in.json
~~~
