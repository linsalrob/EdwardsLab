"""
Create a consensus sequence from a bam file

This is not necessarily trivial, and much of this code comes from
https://github.com/acorg/dark-matter/blob/master/bin/make-consensus.py

You should really use their code.

"""

import os
import sys
import argparse
import logging
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-f', help='reference fasta file', required=True)
    parser.add_argument('-b', help='bam file', required=True)
    parser.add_argument('-o', help='output file')
    parser.add_argument('-v', help='verbose output', action='store_true')

    # completely copied code!
    parser.add_argument(
        '--maskLowCoverage', default=0, type=int,
        help=('Put an N into sites where the coverage is below the specified '
              'cutoff. If you specify a negative numer, masking will be '
              'turned off. Requires --bam.'))

    args = parser.parse_args()

