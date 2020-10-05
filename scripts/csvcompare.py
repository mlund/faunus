#!/usr/bin/env python

import sys

if sys.version_info < (3, 0):
    sys.stdout.write("Sorry, Python 3 og higher required\n")
    sys.exit(1)

import argparse, pandas
from math import fabs

parser = argparse.ArgumentParser(description='Numerically compare two CSV files')
parser.add_argument('--tol', default=0.02, type=float, help='relative error tolerance (default: 0.02)')
parser.add_argument('--small', default=1e-10, type=float, help='always equal if difference is smaller than this (default: 1e-10)')
parser.add_argument('--quiet', '-q', dest='quiet', action='store_true', help='less output')
parser.add_argument('--csv_sep', default=' ', help='CSV separator (default: " ")')
parser.add_argument('file_ref', help='reference file')
parser.add_argument('file_new', help='new file')
args = parser.parse_args()

returncode = 0

def isapprox(a, b):
    return fabs(a - b) < args.small or fabs(a - b) < args.tol * fabs(a)

def isnumber(val):
    ''' filter to compare only ints and floats '''
    return isinstance(val, float) or isinstance(val, int)

def compare(a, b):
    ''' compare ints and floats to relative tolerence '''
    if isnumber(a) and isnumber(b):
        return isapprox(a, b)
    else:
        return a == b

dfs = [pandas.read_csv(f, sep=args.csv_sep, header=None) for f in (args.file_ref, args.file_new)]

if len(dfs) == 2:
    for (ref_col_name, ref_col), (new_col_name, new_col) in zip(*(df.iteritems() for df in dfs)):
        for ref_val, new_val in zip(ref_col, new_col):
            if not compare(ref_val, new_val):
                if returncode == 0:
                    returncode = 1
                if not args.quiet:
                    print('mismatch in col {:2d} {} {}'.format(ref_col_name, ref_val, new_val))
else:
    returncode = 2

sys.exit(returncode)
