"""
Count the sequences in a GFA record
"""

import os
import sys
import argparse
from roblib import bcolors, stream_gfa_sequences
__author__ = 'Rob Edwards'






if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Count sequences in a GFA file ')
    parser.add_argument('-f', help='input file', required=True)
    parser.add_argument('-t', help='tab separated output', action='store_true')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    lens = []
    for i,s in stream_gfa_sequences(args.f):
        if args.v:
            sys.stderr.write(f"{bcolors.GREEN}Parsing {i}{bcolors.ENDC}\n")
        lens.append(len(s))
    
    if len(lens) == 0:
        sys.stderr.write(f"{bcolors.RED}No sequences{bcolors.ENDC}\n")
        sys.exit()

    lens.sort()
    length=sum(lens)
    len_so_far = 0
    n50 = None
    n75 = None

    for i in lens:
        len_so_far += i
        if not n50 and len_so_far >= length * 0.5:
            n50 = i
        if not n75 and len_so_far >= length * 0.75:
            n75 = i

    if args.t:
        print("\t".join(map(str, [args.f, len(lens), length, lens[0], lens[-1], n50, n75])))
    else:
        print("Number of sequences: {}\nTotal length: {}\nShortest: {}\nLongest: {}\nN50: {}\nN75: {}".format(
            len(lens), length, lens[0], lens[-1], n50, n75,
        ))
    
