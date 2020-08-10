"""
Please note: This is a practice for code that will end up in primer-trimming github (that we should rename).

You should use that version
"""

import os
import sys
import argparse
import PyPrinseq

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-f', help='fastq file', required=True)
    parser.add_argument('-l', help='left primers file', required=True)
    args = parser.parse_args()

    PyPrinseq.primertrimming(args.f, args.l, None)
