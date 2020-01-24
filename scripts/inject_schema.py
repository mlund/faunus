#!/usr/bin/env python
#
# This can be used to generate markdown tables
# for topology, moves, energy etc. using the
# information stored in `schema.yml`. It is currently
# inactive as readthedocs would need to call it
#
# The example below can be used to expand the atomlist
# table into `topology.md` by adding the jinja2 variable
# {{ topology_atomlist }}
#
# Should currently be run from the docs/ folder
#
import sys
if sys.version_info < (3, 0):
    sys.stdout.write("Sorry, Python 3 og higher required\n")
    sys.exit(1)

from jinja2 import Template
try:
    try:
        import ruamel.yaml as yaml
    except ImportError:
        import ruamel_yaml as yaml
except ImportError:
    import yaml

def extract_as_lists(schema):
    '''
    Extract keys, types, and description as lists

    Required keywords are suffixed with * and, if
    found, default value is added
    '''
    properties = schema.get("properties", {})
    _keys = []
    _types = []
    _descriptions = []
    if isinstance(properties, dict):
        for key, value in properties.items():
            required = key in schema.get("required", [""])
            _keys.append(key + str("*" if required else ""))
            _descriptions.append(value.get("description", ""))
            _types.append("")
            if "type" in value:
                _types[-1] = value["type"]
                _default = value.get("default", "")
                if _default!="":
                    _keys[-1] += '=' + str(_default).replace(" ", "")
    return zip(_keys, _types, _descriptions)

def print_table(schema, label='Keyword'):
    ''' pretty print schema as markdown table '''
    out = "{:25} | {:7} | {:20}\n".format("Keywords", "Type", "Description")
    out += "{:25} | {:7} | {:40}\n".format(25*"-", 7*"-", 40*"-")
    for k, t, d in extract_as_lists(atomlist):
        out += "{:25} | {:7} | {}\n".format('`'+k+'`', t, d)
    return out

template = Template(sys.stdin.read()) # load template from stdin
with open('schema.yml') as schema_file:
    d = yaml.safe_load(schema_file)["properties"] # read schema file
    atomlist = d["atomlist"]["items"]["additionalProperties"]
    out = template.render(topology_atomlist = print_table(atomlist, '`atomlist`'))
    print(out)
