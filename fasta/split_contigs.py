"""
Split a fasta file of contigs to make them all <= a maximum length (10kb by default).
"""

import os
import sys
import argparse
from roblib import stream_fasta

def split_contigs(inf, outf, length, verbose=False):
    """
    Split the contigs
    :param inf: input fasta file
    :param outf: output fasta file
    :param length: length to split into
    :param verbose: more output
    :return:
    """

    with open(outf, 'w') as out:
        for seqid, seq in stream_fasta(inf, True):
            if verbose:
                sys.stderr.write("{}\n".format(seqid))
            posn=0
            seqcounter=0
            seqidparts = seqid.split(" ")
            while posn < len(seq)-length:
                seqcounter += 1
                if len(seqidparts) == 1:
                    out.write(">{}_{} {}\n".format(seqidparts[0], seqcounter, "".join(seqidparts[1:])))
                else:
                    out.write(">{}_{}\n".format(seqidparts[0], seqcounter))
                out.write("{}\n".format(seq[posn:posn+length]))
                posn += length
            seqcounter += 1
            if len(seqidparts) == 1:
                out.write(">{}_{} {}\n".format(seqidparts[0], seqcounter, "".join(seqidparts[1:])))
            else:
                out.write(">{}_{}\n".format(seqidparts[0], seqcounter))
            out.write("{}\n".format(seq[posn:]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Split contigs to be <= maximum length")
    parser.add_argument('-f', help='fasta file', required=True)
    parser.add_argument('-o', help='output fasta file', required=True)
    parser.add_argument('-l', help='length in bp. Default 10000', type=int, default=10000)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    split_contigs(args.f, args.o, args.l, args.v)