"""
Read a blast output file and count the differences in the sequences.

We have lines that look like this:

    Query  2     CCTGAAATTCTTCGTTCTTCGGACGAGGTGgcagctatcagacagcagagagcgcagcag  61
                 |||||||| || |||||||| ||||||||||| || ||| | ||||||||||| ||||| 
    Sbjct  2     CCTGAAATCCTCCGTTCTTCTGACGAGGTGGCCGCAATCCGCCAGCAGAGAGCACAGCAA  61

So we will look for a line that has Query\s+\d+\s+([A-Z]+)\s+(\d+)\n and a line
that has Sbjct\s+\d+\s+([A-Z]+)\s+(\d+)\n.
Then we will compare the two sequences and count the differences.
"""

import os
import sys
import re
import argparse
from roblib import bcolors
__author__ = 'Rob Edwards'





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='input file', required=True)
    parser.add_argument('-s', help='sort differences', action='store_true')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not os.path.exists(args.f):
        print(f"{bcolors.FAIL}File {args.f} does not exist.{bcolors.ENDC}", file=sys.stderr)
        sys.exit(1)

    if not args.s:
        print(f"{bcolors.WARNING}Warning: Sorting differences is not enabled. Use -s to enable.{bcolors.ENDC}", file=sys.stderr)

    with open(args.f, 'r') as f:
        lines = f.readlines()
    query_pattern = re.compile(r'Query\s+\d+\s+([a-zA-Z\-]+)\s+(\d+)')
    sbjct_pattern = re.compile(r'Sbjct\s+\d+\s+([a-zA-Z\-]+)\s+(\d+)')
    differences = {}

    for i in range(len(lines) - 2):
        query_match = query_pattern.match(lines[i])
        sbjct_match = sbjct_pattern.match(lines[i + 2])
        if query_match and sbjct_match:
            query_seq = query_match.group(1).upper()
            sbjct_seq = sbjct_match.group(1).upper()
            if len(query_seq) != len(sbjct_seq):
                print(f"{bcolors.WARNING}Warning: Sequences are of different lengths.{bcolors.ENDC}")
                continue
            for i in range(len(query_seq)):
                if query_seq[i] != sbjct_seq[i]:
                    if args.s:
                        diff = "".join(sorted((query_seq[i], sbjct_seq[i])))
                    else: 
                        diff = "".join((query_seq[i], sbjct_seq[i]))
                    differences[diff] = differences.get(diff, 0) + 1
    print(f"{bcolors.OKGREEN}Differences found: {len(differences)}{bcolors.ENDC}", file=sys.stderr)

    # print the differences sorted by dictioniary value, highest to lowest
    for d in sorted(differences.items(), key=lambda x: x[1], reverse=True):
        print(f"{bcolors.OKBLUE}{d[0]} : {d[1]} {bcolors.ENDC}")



