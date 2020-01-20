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

if pygments:
    if args.color:
        formatter = Terminal256Formatter
    else:
        formatter = NullFormatter

# Open schema file

def print_table(schema):
    ''' pretty print schema '''
    properties = schema.get("properties", "")
    if (properties!=""):
        print("{:15} | {:8} | {:20}".format("property", "type", "description"))
        print("{:15} | {:8} | {:20}".format(15*"-", 8*"-", 20*"-"))
        for key, value in properties.items():
            required = key in schema.get("required", [""])
            print("{:15} | {:8} | {}".format(key+str("*" if required else ""),
                value["type"], value.get("description", "")))

def validate_input(instance):
    '''
    JSON schema checker

    Will look for schema in faunus/docs/
    '''
    if jsonschema:
        pathname = os.path.dirname(sys.argv[0])
        schemafile = os.path.abspath(pathname+"/../docs") + '/schema.yml'
        if os.path.exists(schemafile):
            with open(schemafile, "r") as f:
                _schema = yaml.safe_load(f)
                error = best_match(Draft7Validator(_schema).iter_errors(instance))
                if error!=None:
                    print( "{}: {}\n".format(error.path, error.message) )
                    print_table(error.schema)
                    sys.exit(1)

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
        validate_input(d)
        out = json.dumps(d, indent=args.indent)
        if pygments:
            out = highlight(out, JsonLexer(), formatter())
        print(out)
    except yaml.parser.ParserError as exception:
        print("input error: invalid json or yaml format", file=sys.stderr)
        print(exception, file=sys.stderr)
