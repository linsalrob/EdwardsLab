"""
Print the locations of the ORFs on the sequence
"""

import os
import sys
import argparse
from roblib import translate_dna, stream_fasta, rc

__author__ = 'Rob Edwards'


def read_sources(sf):
    """
    Read the sources.txt file that has source of ORF call, ORF id
    :param sf: sources.txt file
    :return: dict of ORF/source
    """

    s = {}
    with open(sf, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            s[p[0]] = p[1]
    return s

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Print the locations of the ORFs on the sequence')
    parser.add_argument('-f', help='fasta DNA sequence file', required=True)
    parser.add_argument('-o', help='ORFs file', required=True)
    parser.add_argument('-s', help='sources file that has the source of the ORF calls (default = sources.txt)', default="sources.txt")
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    sources = read_sources(args.s)

    # generate the six frame translations
    seqs = {}
    lengths = {}
    for seqid, seq in stream_fasta(args.f):
        seq = seq.upper()
        lengths[seqid] = len(seq)
        seqs[seqid]["f1"] = translate_dna(seq, args.v)
        seqs[seqid]["f2"] = translate_dna(seq[1:], args.v)
        seqs[seqid]["f3"] = translate_dna(seq[2:], args.v)
        rcseq = rc(seq)
        seqs[seqid]["r1"] = translate_dna(rcseq, args.v)
        seqs[seqid]["r2"] = translate_dna(rcseq[1:], args.v)
        seqs[seqid]["r3"] = translate_dna(rcseq[2:], args.v)

    for orfid, orf in stream_fasta(args.o):
        for s in seqs:
            for fr in ["f1", "f2", "f3", "r1", "r2", "r3"]:
                if orf in seqs[s][fr]:
                    start = seqs[s][fr].index(orf) * 3 # convert from aa to bp!
                    frame = 0
                    end = 0
                    # shit, now we have to figure out frames!
                    if 'f1' == fr:
                        # first frame, this is easy
                        start += 1 # index is 0 indexed!
                        frame = 1
                        end = start + (3 * len(orf)) - 1 # if the aa string includes * this will too.
                    elif 'f2' == fr:
                        start += 2
                        frame = 2
                        end = start + (3 * len(orf)) - 1
                    elif 'f3' == fr:
                        start += 3
                        frame = 3
                        end = start + (3 * len(orf)) - 1
                    elif 'r1' == fr:
                        start = (lengths[s] - start)
                        frame = -1
                        end = start - (3 * len(orf)) + 1
                    elif 'r2' == fr:
                        start = (lengths[s] - (start+1))
                        frame = -2
                        end = start - (3 * len(orf)) + 1
                    elif 'r3' == fr:
                        start = (lengths[s] - (start+2))
                        frame = -3
                        end = start - (3 * len(orf)) + 1


                    print("\t".join(map(str, [sources.get(orfid, 'UNKNOWN'), orfid, s, start, end, frame, lengths[s]])))
                    break
