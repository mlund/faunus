#!/usr/bin/env python

import sys
import json, sys, argparse
try:
    import ruamel_yaml as yaml
except:
    import yaml

try:
    from pygments import highlight
    from pygments.lexers import JsonLexer, YamlLexer
    from pygments.formatters import NullFormatter, Terminal256Formatter
    pygments = True
except:
    pygments = False

parser = argparse.ArgumentParser(description='Convert json to yaml and vice versa')
parser.add_argument('--indent', '-i', dest='indent', default=4, type=int, help='text indentation')
parser.add_argument('--json', '-j', dest='alwaysjson', action='store_true', help='always output to json')
parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
        default=sys.stdin, help='json/yaml input')

if pygments:
    parser.add_argument('--color', '-c', dest='color', action='store_true', help='syntax colored output')

args = parser.parse_args()

if pygments:
    if args.color:
        formatter = Terminal256Formatter
    else:
        formatter = NullFormatter

if sys.version_info < (3, 0):
    sys.stdout.write("Sorry, Python 3 og higher required\n")
    sys.exit(1)

try: # ... to read json
    i = args.infile.read()
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

