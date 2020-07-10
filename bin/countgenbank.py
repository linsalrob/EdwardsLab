"""
Count features in a genbank file or directory of files
"""

import os
import sys
import argparse
from roblib import message, genbank_seqio

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


def count_feats(gbkf, verbose=False):
    if verbose:
        message(f"Reading {gbkf}", "BLUE")

    count = {}
    for seq in genbank_seqio(gbkf):
        for feat in seq.features:
            count[feat.type] = count.get(feat.type, 0) + 1
    return count


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-f', help='genbank file')
    parser.add_argument('-d', help='directory of genbank files')
    parser.add_argument('-t', help='feature type(s) (at least one must be provided)', nargs="+")
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    files = []
    if args.f:
        files.append(args.f)
    if args.d:
        for f in os.listdir(args.d):
            files.append(os.path.join(args.d, f))
    if len(files) == 0:
        message("Fatal. Either -d or -f is required", "RED")

    if len(args.t) == 0:
        message("Fatal. Please provide at least one feature type to count", "RED")

    print("File", end="")
    for t in args.t:
        print(f"\t{t}", end="")
    print()
    for f in files:
        c = count_feats(f, args.v)
        print(f, end="")
        for t in args.t:
            if t in c:
                print(f"\t{c[t]}", end="")
            else:
                print("\t0")
        print()



