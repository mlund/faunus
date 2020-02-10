#!/usr/bin/env python

import sys
if sys.version_info < (3, 0):
    sys.stdout.write("Sorry, Python 3 og higher required\n")
    sys.exit(1)

import os
import json, sys, argparse
import warnings

try:
    jsonschema = True
    from jsonschema import Draft7Validator
    from jsonschema.exceptions import best_match
except ImportError:
    jsonschema = False

try:
    import jinja2
except ImportError:
    warnings.warn("warning: missing jinja2 module")
    jinja2 = None


# API compatibility between pyyaml and ruamel.yaml might break in the future
# https://yaml.readthedocs.io/en/latest/api.html
try:
    try:
        import ruamel.yaml as yaml
    except ImportError:
        # anaconda packs it as a standalone module (underscore instead of dot)
        import ruamel_yaml as yaml
except ImportError:
    import yaml

try:
    pygments = True
    from pygments import highlight
    from pygments.lexers import JsonLexer, YamlLexer
    from pygments.formatters import NullFormatter, Terminal256Formatter
except ImportError:
    pygments = False

parser = argparse.ArgumentParser(description='Convert json to yaml and vice versa')
parser.add_argument('--indent', '-i', dest='indent', default=4, type=int, help='text indentation')
parser.add_argument('--json', '-j', dest='alwaysjson', action='store_true', help='always output to json')
parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='json/yaml input')

if pygments:
    parser.add_argument('--color', '-c', dest='color', action='store_true', help='syntax colored output')

args = parser.parse_args()

red = "\033[91m"
yellow = "\033[93m"
 
if pygments:
    if args.color:
        formatter = Terminal256Formatter
    else:
        formatter = NullFormatter


def human_readable_path(path):
    out=str()
    for i in path:
        if not isinstance(i, int):
            out += str(i)+' -> '
    return out

def eprint(*args, **kwargs):
    ''' print to stderr '''
    print(*args, file=sys.stderr, **kwargs)


def print_table(schema):
    ''' pretty print schema as markdown table '''
    properties = schema.get("properties", "")
    if isinstance(properties, dict):
        eprint(yellow + "Need help, my young apprentice?\n" + red)
        eprint("{:25} | {:7} | {:50}".format(25*"-", 7*"-", 50*"-"))
        eprint("{:25} | {:7} | {:50}".format("Property", "Type", "Description"))
        eprint("{:25} | {:7} | {:50}".format(25*"-", 7*"-", 50*"-"))
        for key, value in properties.items():
            required = key in schema.get("required", [""])
            _property = key + str("*" if required else "")
            if isinstance(value, dict):
                _description = value.get("description", "")
                if "type" in value:
                    _default = value.get("default", "")
                    if _default!="":
                        _property = _property + '=' + str(_default)
                    _type = value["type"]
                    if isinstance(_type, list):
                        _type = '/'.join(_type)
                    eprint("{:25} | {:7} | {}".format('`'+_property+'`', _type, _description))
                else:
                    eprint("{:25} | {:7} | {}".format('`'+_property+'`', 'n/a', _description))

def validate_input(instance):
    ''' JSON schema checker '''
    if jsonschema:
        pathname = os.path.dirname(sys.argv[0]) # location of yason.py
        for subdir in ["/../docs", "/../share/faunus"]: # schema file can be in src or installed
            schemafile = os.path.abspath(pathname+subdir) + '/schema.yml'
            if os.path.exists(schemafile):
                with open(schemafile, "r") as f:
                    _schema = yaml.safe_load(f)
                    error = best_match(Draft7Validator(_schema).iter_errors(instance))
                    if error!=None:
                        eprint( "{}{}\n".format(human_readable_path(error.path), error.message) )
                        print_table(error.schema)
                        sys.exit(1)
                break

try:  # ... to read json
    i = args.infile.read()
    if jinja2:
        # additional files can be used with {% include "file" %}
        dirs = [os.getcwd(), os.path.dirname(os.path.realpath(__file__)) + "/../top"]
        loader = jinja2.FileSystemLoader(dirs)
        env = jinja2.Environment(loader=loader)
        i = env.from_string(i).render()  # render jinja2
        # i = jinja2.Template(i).render() # render jinja2

    d = json.loads(i)
    if "mcloop" in d or "version" in d: 
        validate_input(d)
    if args.alwaysjson:
        if pygments:
            i = highlight(out, JsonLexer(), formatter())
        print(i)
    else:
        out = yaml.safe_dump(d, indent=args.indent, allow_unicode=True)
        if pygments:
            out = highlight(out, YamlLexer(), formatter())
        print(out)
except json.decoder.JSONDecodeError:
    try:  # ... to read yaml
        d = yaml.safe_load(i)  # plain load was deprecated in PyYAML
        if "mcloop" in d or "version" in d: 
            validate_input(d)
        out = json.dumps(d, indent=args.indent)
        if pygments:
            out = highlight(out, JsonLexer(), formatter())
        print(out)
    except yaml.parser.ParserError as exception:
        print("input error: invalid json or yaml format", file=sys.stderr)
        print(exception, file=sys.stderr)
        sys.exit(1)
