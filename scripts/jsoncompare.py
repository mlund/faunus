#!/usr/bin/env python

import sys

if sys.version_info < (3,0):
    sys.stdout.write("Sorry, Python 3 og higher required\n")
    sys.exit(1)

import json, argparse
from math import fabs
try:
    import dictdiffer
except:
    pass

parser = argparse.ArgumentParser(description='Numerically compare two JSON files')
parser.add_argument('--tol', default=0.02, type=float, help='relative error tolerance (default: 0.02)')
parser.add_argument('--small', default=1e-10, type=float, help='always equal if difference is smaller than this (default: 1e-10)')
parser.add_argument('--quiet', '-q', dest='quiet', action='store_true', help='less output')
parser.add_argument('file_ref', help='first file (reference)')
parser.add_argument('file_new', help='second file (new)')
args = parser.parse_args()

returncode = 0

def isapprox(a, b):
    return fabs(a - b) < args.small or fabs(a - b) < args.tol * fabs(a)

def isnumber(key, val):
    ''' filter to compare only ints and floats '''
    if ("time" in key):
        return False
    else:
        return (isinstance(val, float) or isinstance(val, int)) and val != 0

def compare(a, b):
    ''' compare ints and floats in dict to relative tolerence recursively '''
    global returncode
    if isinstance(a, dict):
        for key in set(a.keys()) & set(b.keys()):
            if (key == "groups"): # skip groups
                continue
            if isinstance(a[key], dict):
                compare(a[key], b[key])
            elif isinstance(a[key], list):
                for i, j in zip(a[key], b[key]):
                    compare(i, j)
            elif isnumber(key, a[key]) and isnumber(key, b[key]):
                result = isapprox(a[key], b[key])
                if not result and returncode == 0:
                    returncode = 1
                if not args.quiet:
                    print('{:24} {:>8} {:16.6G} {:16.6G}'.format(key, str(result), a[key], b[key]))

# load two json files to compare
js = [json.load(open(f)) for f in (args.file_ref, args.file_new)]

if len(js) == 2:
    if not args.quiet:
        print('{:24} {:>8} {:>16} {:>16}'.format("keyword", "pass", "reference", "current"))
        print('-'*(24+8+16+16+3))
    compare(*(j for j in js))
    sys.exit(returncode)

sys.exit(1)
