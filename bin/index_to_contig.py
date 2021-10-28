"""
Given a tuple of index1, index2, correlation and a tuple of index, contig

rewrite the correlation to be contig1, contig2, correlation
"""

import os
import sys
import argparse

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-i', help='index file', required=True)
    parser.add_argument('-c', help='correlation file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    i2c = {}
    with open(args.i, 'r') as fin:
        for l in fin:
            if ',' in l:
                p = l.strip().split(',')
            elif "\t" in l:
                p = l.strip().split("\t")
            else:
                sys.stderr.write(f"Neither a comma or tab in {args.i}. What is separator?\n")
                sys.exit(1)

    with open(args.c, 'r') as fin:
        for l in fin:
            if ',' in l:
                p = l.strip().split(',')
            elif "\t" in l:
                p = l.strip().split("\t")
            else:
                sys.stderr.write(f"Neither a comma or tab in {args.c}. What is separator?\n")
                sys.exit(1)
            if p[0] not in i2c:
                sys.stderr.write(f"{p[0]} not found in the index file\n")
                sys.exit(1)
            if p[1] not in i2c:
                sys.stderr.write(f"{p[1]} not found in the index file\n")
                sys.exit(1)
            print(f"{i2c[p[0]]}\t{i2c[p[1]]}\t{p[2]}")