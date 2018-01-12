---
---
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

# Running Simulations

Once compiled, faunus is run from a terminal, taking input either
from `stdin` or a file and can JSON formatted or, via `yason.py` see below,
also YAML formatted:

~~~ bash
faunus < input.json # from json
yason.py minimal.yml | faunus # from yaml
~~~

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

