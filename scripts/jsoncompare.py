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

parser = argparse.ArgumentParser(description='Numerically compare two JSON files')
parser.add_argument('--tol', default=0.02, type=float, help='relative error tolerance (default: 0.02)')
parser.add_argument('--small', default=1e-10, type=float, help='always equal if difference is smaller than this (default: 1e-10)')
parser.add_argument('--quiet', '-q', dest='quiet', action='store_true', help='less output')
parser.add_argument('file1', help='first file (reference)')
parser.add_argument('file2', help='second file (new)')
args = parser.parse_args()

returncode = 0

def equals(a, b):
    if fabs(a-b)<args.small:
        return True
    return fabs(a-b)/fabs(a) < args.tol

def isnumber(key, val):
    ''' filter to compare only ints and floats '''
    if isinstance(val, float) or isinstance(val, int):
        if val!=0:
            if ("time" in key)==False:
                if key!="N_reservoir":
                    return True
    return False

def compare(a, b):
    ''' compare ints and floats in dict to relative tolerence '''
    global returncode
    if isinstance(a, dict):
        for key in set(a.keys()) & set(b.keys()):
            if (key=="groups"): # skip groups
                continue
            if isinstance(a[key], dict):
                compare(a[key], b[key])
            elif isinstance(a[key], list):
                for i, j in zip(a[key], b[key]):
                    compare(i,j)
            elif isnumber(key, a[key]) and isnumber(key, b[key]):
                result = equals(a[key], b[key])
                if result==False and returncode==0:
                    returncode = 1
                if not args.quiet:
                    print('{:24} {:>8} {:16.6G} {:16.6G}'.format(key, str(result), a[key], b[key]))

# load two json files to compare
d = [json.load(open(f)) for f in [args.file1, args.file2]]

if (len(d)==2):
    if args.quiet==False:
        print('{:24} {:>8} {:>16} {:>16}'.format("keyword", "pass", "reference", "current"))
        print('-'*(24+8+16+16+3))
    compare(*[i for i in d])
    sys.exit(returncode)

sys.exit(1)

