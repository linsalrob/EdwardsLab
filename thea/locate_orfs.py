"""
Print the locations of the ORFs on the sequence
"""

import os
import sys
import argparse
from roblib import translate_dna, stream_fasta, rc

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Print the locations of the ORFs on the sequence')
    parser.add_argument('-f', help='fasta DNA sequence file', required=True)
    parser.add_argument('-o', help='ORFs file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    # generate the six frame translations
    seqs = {}
    for seqid, seq in stream_fasta(args.f):
        if args.v:
            sys.stderr.write("Translating frame 1\n")
        seqs["{}_f1".format(seqid)] = translate_dna(seq, args.v)
        if args.v:
            sys.stderr.write("Translating frame 2\n")
        seqs["{}_f2".format(seqid)] = translate_dna(seq[1:], args.v)
        if args.v:
            sys.stderr.write("Translating frame 3\n")
        seqs["{}_f3".format(seqid)] = translate_dna(seq[2:], args.v)
        rcseq = rc(seq)
        if args.v:
            sys.stderr.write("Translating frame -1\n")
        seqs["{}_r1".format(seqid)] = translate_dna(rcseq, args.v)
        if args.v:
            sys.stderr.write("Translating frame -2\n")
        seqs["{}_r2".format(seqid)] = translate_dna(rcseq[1:], args.v)
        if args.v:
            sys.stderr.write("Translating frame -3\n")
        seqs["{}_r3".format(seqid)] = translate_dna(rcseq[2:], args.v)

    for orfid, orf in stream_fasta(args.o):
        for s in seqs:
            if orf in seqs[s]:
                start = seqs[s].index(orf)
                frame = 0
                # shit, now we have to figure out frames!
                if s.endswith('_f1'):
                    # first frame, this is easy
                    start += 1 # index is 0 indexed!
                    frame = 1
                    end = start + (3 * len(orf))
                elif s.endswith('_f2'):
                    start += 2
                    frame = 2
                    end = start + (3 * len(orf))
                elif s.endswith('_f3'):
                    start += 2
                    frame = 3
                    end = start + (3 * len(orf))
                elif s.endswith('_r1'):
                    start = (len(seqs[s]) - start) + 1
                    frame = -1
                    end = start - (3 * len(orf))
                elif s.endswith('_r2'):
                    start = (len(seqs[s]) - start) + 2
                    frame = -2
                    end = start - (3 * len(orf))
                elif s.endswith('_r3'):
                    start = (len(seqs[s]) - start) + 3
                    frame = -3
                    end = start - (3 * len(orf))

                print("\t".join(map(str, [orfid, start, end, frame])))
                break