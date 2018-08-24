"""
Third attempt to read a set of sequences and a list of prophages and make sequences without the prophages.

This is designed for the output from the revised phispy
"""

import os
import sys
import argparse
from operator import itemgetter
from roblib import stream_fasta

def read_phage_locations(locf, verbose=False):
    """
    Read the phage locations. Expects tab separated text of contig, start, stop, number of genes
    :param locf: the locations file
    :param verbose: more output
    :return: a dict of all the contigs that have prophages, and then tuples of [start, stop, number of genes]
    """

    contigs = {}
    with open(locf, 'r') as f:
        for l in f:
            if p.endswith('Status'):
                continue
            p=l.strip().split("\t")
            if p[1] > p[2]:
                (p[1], p[2]) = (p[2], p[1])
            if p[0] not in contigs:
                contigs[p[0]] = []
            contigs[p[0]].append((int(p[1]), int(p[2]), p[3]))
    return contigs


def read_sequence(conf, verbose=False):
    """
    Read the contigs file for this genome and return it
    :param conf: the contigs file
    :param verbose:
    :return: a dict of contig/sequence
    """

    seqs = {}
    for seqid, seq in stream_fasta(conf, whole_id=False):
        seqs[seqid] = seq
    return seqs

def write_seqs(seqs, locations, phagef, nonphagef, verbose=False):
    """
    Write the sequences out
    :param phagef: file to write the prophages to
    :param nonphagef: file to write the nonprophages to
    :param seqs: DNA sequences
    :param locations: locations of phage regions
    :param verbose: more output
    :return:
    """

    with open(phagef, 'w') as phageout:
        with open(nonphagef, 'w') as nonphageout:
            for s in seqs:
                if s in locations:
                    ses = sorted(locations[s], key=itemgetter(0))
                    posn = 0
                    for start, end, numgenes in ses:
                        nonphagef.write(">{}\n{}\n".format("_".join(map(str, [s, posn, start-1])),
                                                           seqs[s][posn:start]))
                        phagef.write(">{} [predicted by PhiSpy to have {} genes]\n{}\n".
                                     format("_".join(map(str, [s, start, end])), numgenes, seqs[s][start,end]))
                        posn=end+1
                    nonphagef.write(">{}\n{}\n".format("_".join(map(str, [s, posn, len(seqs[s])])),
                                                       seqs[s][posn:]))
                else:
                    nonphagef.write(">{}\n{}\n".format("_".join(map(str, [s, 0, len(seqs[s])])),
                                                       seqs[s]))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Separate prophage and non prophage DNA")
    parser.add_argument('-l', help='locations in tsv format directory', required=True)
    parser.add_argument('-c', help='contigs directory', required=True)
    parser.add_argument('-n', help='name of directory to write non-prophage sequences', required=True)
    parser.add_argument('-p', help='name of directory to write prophage sequences', required=True)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    for d in ['prophage_seqs', 'nonprophage_seqs']:
        if os.path.exists(d):
            sys.stderr.write("ERROR: {} exists. Not overwriting\n".format(d))
            sys.exit(-1)
    os.mkdir('prophage_seqs')
    os.mkdir('nonprophage_seqs')

    for f in os.listdir(args.l):
        if os.path.exists(os.path.join(args.c. f)):
            loc = read_phage_locations(os.path.join(args.l. f))
            dna = read_sequence(os.path.join(args.c. f))
            write_seqs(dna, loc, os.path.join("prophage_seqs", f), os.path.join("nonprophage_seqs", f))