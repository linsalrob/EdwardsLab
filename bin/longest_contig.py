"""
Print the length of the longest contig for each file in a directory of fasta files.
"""

import os
import sys
import argparse
from roblib import read_fasta, bcolors

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Print the length of the longest contig for each file in a directory of fasta files')
    parser.add_argument('-d', help='Directory of fasta files', required=True)
    parser.add_argument('-f', help='fasta file to write the longest sequence to')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    endings = {'.fna', '.fasta', '.fa'}

    for f in os.listdir(args.d):
        longest = [0, None, None]
        isfasta = False
        for e in endings:
            if f.endswith(e):
                isfasta = True
                break
        if not isfasta:
            if args.v:
                sys.stderr.write(f"{bcolors.PINK}Don't think {f} is a fasta file. Skipped\n{bcolors.ENDC}")
            continue
        if args.v:
            sys.stderr.write(f"{bcolors.GREEN}{f}{bcolors.ENDC}\n")
        fa = read_fasta(os.path.join(args.d, f))
        for x in fa:
            if len(fa[x]) > longest[0]:
                longest = [len(fa[x]), x, fa[x]]
        if 0 == longest[0]:
            continue
        print("{}\t{}".format(f, longest[0]))
        if args.f:
            with open(args.f, 'a') as out:
                out.write(f">{longest[1]} [from {f}]\n{longest[2]}\n")



