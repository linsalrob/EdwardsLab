#!/usr/bin/env python

"""
Renumber multiple fasta files, starting with either number 1 or a number provided on the command line.
At the end, we output the last number written.

Note that this uses the roblib library that is available from this github repository.

This version will combine multiple fasta files into a single output. 
"""

import os
import sys
import argparse
from roblib import read_fasta, bcolors

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='input file', nargs="*")
    parser.add_argument('-d', help='directory of input files', nargs="*")
    parser.add_argument('-o', help='output file', required=True)
    parser.add_argument('-i', help='id mapping file', required=True)
    parser.add_argument('-n', help='number to start from. Default: 1', type=int, default=1)
    parser.add_argument('-x', help='maximum number of sequences to write out', type=int)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()


    counter = args.n - 1

    if not args.f and not args.d:
        sys.stderr.write(f"{bcolors.RED}FATAL: Please supply either -d or -f options{bcolors.ENDC}\n")
        sys.exit(-1)
    
    idmap = open(args.i, 'w')
    out = open(args.o, 'w')

    for f in args.f:
        if args.v:
            sys.stderr.write(f"{bcolors.GREEN}Reading {f}{bcolors.ENDC}\n")
        fa = read_fasta(f)
        for id in fa:
            counter += 1
            out.write(">{}\n{}\n".format(counter, fa[id]))
            idmap.write("{}\t{}\t{}".format(f, id, counter))
            if args.x and (counter - (args.n-2)) > args.x:
                break

    if args.d:
        for d in args.d:
            if args.v:
                sys.stderr.write(f"{bcolors.GREEN}Reading {d}{bcolors.ENDC}\n")
            for f in os.listdir(d):
                if args.v:
                    sys.stderr.write(f"{bcolors.BLUE}\tReading {f}{bcolors.ENDC}\n")
                fa = read_fasta(os.path.join(d, f))
                for id in fa:
                    counter += 1
                    out.write(">{}\n{}\n".format(counter, fa[id]))
                    idmap.write("{}\t{}\t{}".format(f, id, counter))
                    if args.x and (counter - (args.n - 2)) > args.x:
                        break


    idmap.close()
    out.close()
    print("The last ID written to the file {} was {}".format(args.o, counter))

