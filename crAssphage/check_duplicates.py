"""
Parse the NCBI file and check for duplicates. We ignore sample name, title, bioproject accession and description
"""

import os
import sys
import argparse
from roblib import bcolors
__author__ = 'Rob Edwards'


def parse_file(fi, verbose):
    """
    Parse a file
    :param: fi: File to parse
    :param: verbose: more output
    """

    ignore = set()
    seen = {}
    with open(fi, 'r') as f:
        for l in f:
            p=l.strip().split("\t")
            if not ignore:
                for c in ['sample_name', 'title', 'bioproject_accession', 'description']:
                    if c in p:
                        ignore.add(p.index(c))
                continue
            s=""
            t=""
            for i in (range(len(p))):
                if i in ignore:
                    t += f"{bcolors.BLUE}{p[i]}{bcolors.WHITE}"+" | "
                    continue
                t+=p[i]+" | "
                s+=p[i]+" | "
            if s in seen:
                sys.stderr.write(f"{bcolors.RED}DUPLICATE:{bcolors.ENDC}: {seen[s]}\n")
            seen[s]=t



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='input file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    parse_file(args.f, args.v)


