"""
Convert a triple of [name, lat, lon] to a matrix of pairwise distances
"""

import os
import sys
import argparse
from roblib import geography

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert a triple of [name, lat, lon] to a matrix of pairwise distances')
    parser.add_argument('-f', help='input file with [name, lat, lon]', required=True)
    parser.add_argument('-i', help='ignore the first line', action='store_true')
    args = parser.parse_args()

    first = args.i
    loc = {}
    with open(args.f, 'r') as f:
        for l in f:
            if first:
                first = False
                continue
            p = l.strip().split("\t")
            loc[p[0]]=[float(p[1]), float(p[2])]

    names = list(loc.keys())
    names.sort()
    sys.stdout.write("\t" + "\t".join(names))
    sys.stdout.write("\n")
    for i in range(len(names)):
        fn = names[i]
        sys.stdout.write(fn)
        for j in range(len(loc)):
            tn = names[j]
            if fn == tn:
                sys.stdout.write("\t0")
                continue
            #sys.stdout.write("{}, {} :: {}, {}".format(loc[fn][0], loc[fn][1], loc[tn][0], loc[tn][1]))
            dist = geography.latlon2distance(loc[fn][0], loc[fn][1], loc[tn][0], loc[tn][1])
            sys.stdout.write("\t{}".format(dist))
            #sys.stdout.write("\n")
        sys.stdout.write("\n")

