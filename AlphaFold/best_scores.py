"""
Read all the ranking_debug.json files and report the best model and its pLDDT score
"""

import os
import sys
import argparse
import json

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-d', help='directory with ranking_debug.json', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    rdj = os.path.join(args.d, "ranking_debug.json")
    if not os.path.exists(rdj):
        sys.stderr.write(f"Error: {rdj} not found\n")
        sys.exit(1)

    f = open(rdj, 'r')
    d = json.load(f)
    bm = d['order'][0]
    pl = d['plddts'][bm]
    print("\t".join(map(str, [args.d, bm, pl])))