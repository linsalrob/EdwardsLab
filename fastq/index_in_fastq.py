"""
Print the location of a substring in a fastq file
"""

import os
import sys
import argparse
from roblib import stream_fastq

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

def print_locations(fastqf, s, pl):
    """
    Print the location of s in all reads in fastqf
    :param fastqf:
    :param s:
    :param pl: print the sequence length
    :return:
    """

    for seqid, header, seq, qual in stream_fastq(fastqf):
        r = seq.find(s)
        while r > -1:
            if pl:
                print(f"{seqid}\t{r}\t{len(seq)}")
            else:
                print(f"{seqid}\t{r}")
            r += 1
            r = seq.find(s, r)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-f', help='fastq file', required=True)
    parser.add_argument('-s', help='substring to look for', required=True)
    parser.add_argument('-l', help='print sequence length', action='store_true')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    print_locations(args.f, args.s, args.l)