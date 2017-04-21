#!/usr/bin/env python

import sys
import json, sys, argparse
#from ruamel import yaml  # recommended, but somehow doesn't work w. anaconda
import yaml

parser = argparse.ArgumentParser(description='Convert json to yaml and vice versa')
parser.add_argument('--indent', '-i', dest='indent', default=4, type=int, help='text indentation')
parser.add_argument('--json', '-j', dest='alwaysjson', action='store_true', help='always output to json')
parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
        default=sys.stdin, help='json/yaml input')
args = parser.parse_args()

if sys.version_info < (3, 0):
    sys.stdout.write("Sorry, Python 3.x required\n")
    sys.exit(1)

try: # ... to read json
    i = args.infile.read()
    d = json.loads( i )
    if args.alwaysjson:
        print( i )
    else:
        print( yaml.safe_dump(d, indent=args.indent, allow_unicode=True ) )
except:
    try: # ... to read yaml
        d = yaml.load( i )
        print( json.dumps(d, indent=args.indent) )
    except:
        print("input error: invalid json or yaml format")

