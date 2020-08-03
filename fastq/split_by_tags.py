"""
Split fastq files by sequence tags. You need to provide the tags in this code!
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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-f', help='fastq file', required=True)
    parser.add_argument('-o', help='output directory', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    tags = ['GGTTCACTTGAGACAC','CTTGAGACAC']

    os.makedirs(args.o, exist_ok=True)
    outfiles = {'none': open(os.path.join(args.o, "none.fastq"), 'w')}

    for seqid, header, seq, qual in stream_fastq(args.f):
        written = False
        for t in tags:
            if t in seq and seq.index(t) < 25:
                tag = seq[0:seq.index(t)+len(t)]
                if tag not in outfiles:
                    outfiles[tag] = open(os.path.join(args.o, f"{tag}.fastq"), 'w')
                outfiles[tag].write(f"@{header}\n{seq}\n+\n{qual}\n")
                written = True
                break
        if not written:
            outfiles['none'].write(f"@{header}\n{seq}\n+\n{qual}\n")

    for t in outfiles:
        outfiles[t].close()
