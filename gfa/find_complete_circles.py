"""
Parse the gfa graph and find all the complete circles in the graph.
"""

import os
import sys
import argparse
import gzip
from roblib import bcolors
__author__ = 'Rob Edwards'


def find_complete_circles(gfa, verbose=False):
    """
    Read the gfa file and find all the complete circles in the graph.
    A circle is defined as a path that starts and ends at the same node, and at the moment
    we only consider circles that are single contigs
    """
    
    opener=open
    if gfa.endswith('.gz'):
        opener=gzip.open


    fwd = set()
    rev = set()
    seqlen = {}
    connections = {}
    with opener(gfa, 'rt') as f:
        for l in f:
            if l.startswith('S'):
                p = l.strip().split('\t')
                seqlen[p[1]] = len(p[2])

            if l.startswith('L'):
                p = l.strip().split('\t')
                if len(p) != 7:
                    print(f"Error: line {l} does not have 7 fields")
                    continue
                if p[1] not in connections:
                    connections[p[1]] = set()
                if p[3] not in connections:
                    connections[p[3]] = set()
                connections[p[1]].add(p[3])
                connections[p[3]].add(p[1])
                if p[1] == p[3]:
                    if p[2] == '+' and p[4] == '+':
                        fwd.add(p[1])
                    elif p[2] == '-' and p[4] == '-':
                        rev.add(p[1])
                    else:
                        # print(f"{l[1]} has +/- orientation", file=sys.stderr)
                        continue

    for i in fwd:
        if len(connections[i]) == 1:
            print(f"{bcolors.OKGREEN}Circle: {i} {connections[i]} (dir +) (len: {seqlen[i]}){bcolors.ENDC}")
    for i in rev:
        if len(connections[i]) == 1:
            print(f"{bcolors.BLUE}Circle: {i} {connections[i]} (dir -) (len: {seqlen[i]}){bcolors.ENDC}")


    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-g', '--gfa', help='gfa input file', required=True)
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not os.path.exists(args.gfa):
        print(f"{bcolors.FAIL}Error: {args.gfa} does not exist{bcolors.ENDC}")
        sys.exit(1)

    find_complete_circles(args.gfa, args.verbose)


