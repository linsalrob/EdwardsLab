"""
Parse a cd-hit output file and a fasta file and make a directory of clusters.

NOTE: In this case the sequence id's are ALL numbers (i.e. we have used renumber_fasta.py on them)
"""

import os
import sys
import argparse
from roblib import sequences
import re

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse a cd-hit output file and a fasta file and make a directory of clusters')
    parser.add_argument('-f', help='fasta file for DNA sequences', required=True)
    parser.add_argument('-c', help='cd-hit output file with the clusters', required=True)
    parser.add_argument('-o', help='output directory', required=True)
    args = parser.parse_args()

    fa = sequences.read_fasta(args.f)

    if os.path.exists(args.o):
        sys.stderr.write("ERROR: OUTPUT DIRECTORY EXISTS")
        sys.exit(-1)

    os.mkdir(args.o)

    cluster = None
    out = None
    with open(args.c, 'r') as fin:
        for l in fin:
            if l.startswith('>Cluster'):
                m=re.match('>Cluster\s+(\d+)', l)
                cluster = m.group(1)
                if out:
                    out.close()
                out = open(os.path.join(args.o, "{}.fa".format(cluster)), 'w')
            else:
                m=re.search('>(\d+)\.\.\.', l)
                out.write(">{}\n{}".format(m.groups(1), fa[m.groups(1)]))

    out.close()
