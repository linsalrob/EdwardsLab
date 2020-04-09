"""
Extract the ends of sequences to a separate file
"""

import os
import sys
import argparse
from roblib import colours, stream_fastq

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, BacterialProteins'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


def split_fastq(fqf, outdir, frac, verbose=False):
    """
    Split a fastq file
    :param fqf: fastq file
    :param outdir: output directory to write all the files to
    :param frac: fraction of the sequence for each end
    :param verbose: more output
    :return: nothing
    """

    if not os.path.exists(outdir):
        os.path.mkdir(outdir)

    for seqid, header, seq, qual in stream_fastq(fqf):
        with open (os.path.join(outdir, seq + ".left.fna"), 'w') as out:





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-q', help='input fastq file', required=True)
    parser.add_argument('-x', help='fraction of the read from each end to put into separate files (<=1)', type=float)
    parser.add_argument('-o', help='output file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()




