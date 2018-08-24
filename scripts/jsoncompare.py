#!/usr/bin/env python

import sys

if sys.version_info<(3,0):
    sys.stdout.write("Sorry, Python 3 og higher required\n")
    sys.exit(1)

import json, argparse
from math import fabs
try:
    import dictdiffer
except:
    pass

parser = argparse.ArgumentParser(description='Nulerically compare two JSON files')
parser.add_argument('--tol', default=0.05, type=float, help='relative error tolerance (default: 0.05)')
parser.add_argument('--quiet', '-q', dest='quiet', action='store_true', help='less output')
parser.add_argument('file1', help='first file')
parser.add_argument('file2', help='second file')
args = parser.parse_args()

def equals(a, b):
    return fabs(a-b)/a < args.tol

def isnumber(key, val):
    ''' filter to compare only ints and floats '''
    if isinstance(val, float) or isinstance(val, int):
        if val!=0:
            if ("time" in key)==False:
                return True
    return False

def compare(a, b):
    ''' compare ints and floats in dict '''
    if isinstance(a, dict):
        for key in set(a.keys()) & set(b.keys()):
            if isinstance(a[key], dict):
                compare(a[key], b[key])
            elif isinstance(a[key], list):
                for i, j in zip(a[key], b[key]):
                    compare(i,j)
            elif isnumber(key, a[key]) and isnumber(key, b[key]):
                if args.quiet:
                    return equals(a[key], b[key])
                else: print('{:20} {:10} {:10} {}'.format(key, a[key], b[key], equals(a[key], b[key])))
    return False

# load two json files to compare
d = [json.load(open(f)) for f in [args.file1, args.file2]]

if (len(d)==2):
    if args.quiet==False:
        print('{:20} {:>10} {:>10} {}'.format("keyword", "file1", "file2", "pass"))
    compare(*[i for i in d])
    sys.exit(0)

sys.exit(1)

