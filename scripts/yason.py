#!/usr/bin/env python

import sys
if sys.version_info < (3, 0):
    sys.stdout.write("Sorry, Python 3 og higher required\n")
    sys.exit(1)

import os
import json, sys, argparse
import warnings

try:
    import jinja2
except ImportError:
    warnings.warn("warning: missing jinja2 module");
    jinja2 = None 

try:
    import ruamel_yaml as yaml
    warnings.simplefilter('ignore', yaml.error.UnsafeLoaderWarning)
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

try: # ... to read json
    i = args.infile.read()
    if jinja2:
        # additional files can be used with {% include "file" %}
        dirs = [os.getcwd(), os.path.dirname(os.path.realpath(__file__)) + "/../top"]
        loader = jinja2.FileSystemLoader(dirs)
        env = jinja2.Environment(loader=loader)
        i = env.from_string(i).render() # render jinja2
        #i = jinja2.Template(i).render() # render jinja2 

    d = json.loads( i )
    if args.alwaysjson:
        if pygments:
            i = highlight( out, JsonLexer(), formatter() )
        print( i )
    else:
        out = yaml.safe_dump(d, indent=args.indent, allow_unicode=True )
        if pygments:
            out = highlight( out, YamlLexer(), formatter() )
        print( out )
except:
    try: # ... to read yaml
        d = yaml.load( i )
        out = json.dumps(d, indent=args.indent)
        if pygments:
            out = highlight(out, JsonLexer(), formatter() )
        print(out)
    except:
        print("input error: invalid json or yaml format")

