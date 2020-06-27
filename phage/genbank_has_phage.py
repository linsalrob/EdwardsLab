"""
Test a directory of genbank files and note whether they have the is_phage qualifier for their genomes
"""

import os
import sys
import argparse

from roblib import genbank_seqio, message

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-d', help='directory of genbank files', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    for f in os.listdir(args.d):
        if args.v:
            message(f"Reading {f}", "GREEN")
        pc = 0
        for s in genbank_seqio(os.path.join(args.d, f)):
            for feat in s.features:
                if 'is_phage' in feat.qualifiers:
                    pc += 1
        print(f"{f}\t{pc}")