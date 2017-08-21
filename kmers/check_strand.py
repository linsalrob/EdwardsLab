"""
Given a fasta file of sequences check to see which sequences should be reverse complemented to match best to the first
sequence
"""

import os
import sys
import argparse
from roblib import sequences
from roblib import dna
from scipy.spatial.distance import cdist
import numpy as np

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Given a fasta file of sequences check to see which sequences should be reverse complemented to match best to the first sequence')
    parser.add_argument('-f', help='fasta file of sequences', required=True)
    parser.add_argument('-k', help="k-mer size (default=7)", default=7)
    args = parser.parse_args()

    # for the first sequence we need to store its information
    first = {}
    donefirst = False

    for (seqid, seq) in sequences.stream_fasta(args.f):
        counts = dna.kmers(seq, args.k)
        if not donefirst:
            first['fwd'] = counts

        keys = set(counts.keys())
        keys.update(first['fwd'].keys())
        keys = list(keys)

        std = [first['fwd'].get(k, 0) for k in keys]
        tst = [counts.get(k,0) for k in keys]

        fwddist = cdist(np.array(std), np.array(tst), metric="euclidean")


        counts = dna.kmers(dna.rc(seq), args.k)
        if not donefirst:
            first['rev'] = counts

        keys = set(counts.keys())
        keys.update(first['rev'].keys())
        keys = list(keys)

        std = [first['rev'].get(k, 0) for k in keys]
        tst = [counts.get(k,0) for k in keys]

        revdist = cdist(np.array(std), np.array(tst), metric="euclidean")

        rc = "no"
        if revdist < fwddist:
            rc = "yes"

        print("{}\t{}\t{}\t{}".format(seqid, fwddist, revdist, rc))