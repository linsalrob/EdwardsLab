"""

We're scoring the alignments by looking at one change, pairs of changes, triples of changes, etc.

"""

import os
import sys
import argparse
from roblib import read_fasta

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='fasta alignment file', required=True)
    parser.add_argument('-m', help='minimum size of (mis)matches to count (default=1)', default=1, type=int)
    parser.add_argument('-s', help='maximum size of (mis)matches to count (default=3)', default=3, type=int)
    args = parser.parse_args()

    fa = read_fasta(args.f)

    for sid in fa:
        fa[sid] = fa[sid].upper()

    # score all pairwise comparisons
    seqs = list(fa.keys())
    for i, seqid1 in enumerate(seqs):
        for j in range(i+1, len(seqs)):
            seqid2 = seqs[j]
            s1 = fa[seqid1]
            s2 = fa[seqid2]
            scores = []
            for sublen in range(args.m, args.s+1):
                sim = 0
                diff = 0
                posn = 0
                while posn + sublen < len(s1):
                    s1substr = s1[posn:posn+sublen]
                    s2substr = s2[posn:posn+sublen]

                    if s1substr == s2substr:
                        sim += 1
                        continue
                    alldiff = True
                    for substr_i in range(len(s1substr)):
                        if s1substr[substr_i] == s2substr[substr_i]:
                            alldiff = False
                    if alldiff:
                        diff += 1
                    posn += sublen
                # scores.append(1.0 * diff / (sim + diff))
                scores.append(diff)
            sys.stdout.write("{}\t{}\t".format(seqid1, seqid2))
            sys.stdout.write("\t".join(map(str, scores)))
            sys.stdout.write("\n")

