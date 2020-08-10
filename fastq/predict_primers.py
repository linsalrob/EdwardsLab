"""
Please note: This is a practice for code that will end up in primer-trimming github (that we should rename).

You should use that version

"""

import os
import sys
import argparse
import faulthandler
import PyPrinseq

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-f', help='file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    kmerlen = 10
    minpercent = 10.0
    fasta_ouput = False
    three_prime = False
    print_kmer_counts = False
    print_abundance = False
    print_short = False
    debug = False

    primers = PyPrinseq.primerpredict(args.f, 8, 1, three_prime, debug)

    print(f"There are {len(primers)} primers")

    for i,p in enumerate(primers):
        print(f"Python primer {i} :  {p}")
