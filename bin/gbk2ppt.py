"""
Convert a genbank file to a ppt file
"""

import os
import sys
import argparse
from genbank_to_ppt import convert_genbank


__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert a genbank file to a ppt file")
    parser.add_argument('-f', help='genbank file to convert', required=True)
    parser.add_argument('-o', help='output file')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if args.o:
        r = convert_genbank(args.f, False, args.v)
        with open(args.o, 'w') as out:
            for l in r:
                out.write("\t".join(l))
                out.write("\n")
    else:
        convert_genbank(args.f, True, args.v)
