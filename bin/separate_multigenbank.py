"""
Separate a single genbank file with lots of entries into a directory of single files.
"""

import os
import sys
import argparse
import gzip
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Separate a single genbank file with lots of entries into a directory of single files")
    parser.add_argument('-f', help='genbank file to separate', required=True)
    parser.add_argument('-d', help='output directory to write to', required=True)
    parser.add_argument('-v', help='verbose output')
    args = parser.parse_args()

    if not os.path.exists(args.d):
        os.makedirs(args.d)

    if args.f.endswith('.gz'):
        fin = gzip.open(args.f, 'rt')
    else:
        fin = open(args.f, 'r')

    out = None

    for l in fin:
        if '//' == l:
            if out:
                out.write(l)
                out.close()
                out = None
            continue
        if 'LOCUS' in l:
            if out:
                out.close()
            m=re.match('LOCUS\s+(\S+)', l)
            outfilename = m.groups()[0]
            if not outfilename:
                sys.stderr.write("FATAL: Could not parse a filename from {}".format(l))
                sys.exit(-1)
            outfilename += ".gbk"
            if args.v:
                sys.stderr.write("Writing to {}\n".format(outfilename))
            out = open(os.path.join(args.d, outfilename), 'w')
        if out:
            out.write(l)
        elif args.v:
            sys.stderr.write("SKIPPED: {}".format(l))
