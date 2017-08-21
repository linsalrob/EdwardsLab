"""
Given a fasta file of sequences check to see which sequences should be reverse complemented to match best to the first
sequence
"""

import os
import sys
import argparse
from roblib import sequences
from roblib import dna
import numpy as np

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Given a fasta file of sequences check to see which sequences should be reverse complemented to match best to the first sequence')
    parser.add_argument('-f', help='fasta file of sequences', required=True)
    parser.add_argument('-k', help="k-mer size (default=7)", default=7, type=int)
    parser.add_argument('-d', help='debug output', action="store_true")
    args = parser.parse_args()

    # for the first sequence we need to store its information
    first = {}
    donefirst = False

    for (seqid, seq) in sequences.stream_fasta(args.f):
        fwdcounts = dna.kmers(seq, args.k)
        if not donefirst:
            first['fwd'] = fwdcounts

        keys = set(fwdcounts.keys())
        keys.update(first['fwd'].keys())

        revcounts = dna.kmers(dna.rc(seq), args.k)
        if not donefirst:
            first['rev'] = revcounts

        keys.update(set(revcounts.keys()))
        keys.update(first['rev'].keys())
        keys = list(keys)

        if args.d and not donefirst:
            for k in keys:
                print("{}\t{}\t{}".format(k, first['fwd'].get(k, 0), first['rev'].get(k, 0)))

        fwdstd = [first['fwd'].get(k, 0) for k in keys]
        fwdtst = [fwdcounts.get(k, 0) for k in keys]
        # revstd = [first['rev'].get(k, 0) for k in keys]
        revtst = [revcounts.get(k, 0) for k in keys]

        if args.d and donefirst:
            print(" ".join(map(str, fwdstd[0:60])))
            print(" ".join(map(str, fwdtst[0:60])))
            #print(" ".join(map(str, revstd[0:60])))
            print(" ".join(map(str, revtst[0:60])))

        fwddist = np.linalg.norm(np.array(fwdstd) - np.array(fwdtst))
        # revdist = np.linalg.norm(np.array(revstd) - np.array(revtst))
        revdist = np.linalg.norm(np.array(fwdstd) - np.array(revtst))

        rc = "no"
        if revdist < fwddist:
            rc = "yes"

        print("{}\t{}\t{}\t{}".format(seqid, fwddist, revdist, rc))

        if args.d and donefirst:
            sys.exit(0)
        
        donefirst = True

