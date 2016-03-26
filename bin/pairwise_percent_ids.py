"""
Calculate the pairwise percent identities from protein alignments. We only consider identical sequences
in the alignment
"""

import os
import sys

import argparse


def read_fasta(filename):
    seq = {}
    seqid = None
    prot = ''
    with open(filename, 'r') as f:
        for l in f:
            if l.startswith('>'):
                if seqid:
                    seq[seqid] = prot
                prot = ''
                seqid = l.strip().replace('>', '')
            else:
                prot += l.strip()
    seq[seqid] = prot
    return seq


def pairwise(seqs):
    allseqs = seqs.keys()
    allseqs.sort()
    for i in range(len(allseqs)):
        for j in range(len(allseqs)):
            maxp = max(len(seqs[allseqs[i]]), len(seqs[allseqs[j]]))
            same = 0
            diff = 0
            for p in range(maxp):
                if seqs[allseqs[i]][p] == '-' and seqs[allseqs[j]][p] == '-':
                    pass
                if seqs[allseqs[i]][p] == seqs[allseqs[j]][p]:
                    same += 1
                else:
                    diff += 1
            print("{}\t{}\t{}".format(allseqs[i], allseqs[j], (1.0 * same/(same+diff))*100))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate the percent pairwise ID between all pairs of sequences")
    parser.add_argument('-f', help='protein fasta file', required=True)
    args = parser.parse_args()

    sq = read_fasta(args.f)
    pairwise(sq)
