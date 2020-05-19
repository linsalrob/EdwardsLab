"""
Count the leading k-mers in a fastq file. This is useful to see if the 
k-mers should be trimmed
"""

import os
import sys
import argparse
from roblib import bcolors,stream_fastq
__author__ = 'Rob Edwards'

k = {"A" : 0, "G" : 1, "C" : 2, "T" : 3,
 "a" : 0, "g" : 1, "c" : 2, "t" : 3}

def count_kmers(fqf, kmer, verbose=False):
    """ 
    Count hte frequency of bases in the first k-mer bp of a sequences
    :param fqf: fastq file
    :param kmer: length to count
    :param verbose: more output
    """

    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Reading {fqf}{bcolors.ENDC}\n");
    counts = [[0,0,0,0] for x in range(kmer)]
    for sid, seqid, seq, qual in stream_fastq(fqf):
        if 'N' in seq  or 'n' in seq:
            continue
        seq = seq[::-1] # reverse the sequence!
        for x in range(kmer):
            try:
                counts[x][k[seq[x]]] += 1
            except KeyError as e:
                base = e.args[0]
                if base.upper() != "N":
                    sys.stderr.write(f'{bcolors.PINK}Unknown base {base}{bcolors.ENDC}\n')
    counts.reverse()
    return counts

def predict(counts, cutoff=0.5, verbose=False):
    """
    What sequence should be removed?
    """
    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Predicting{bcolors.ENDC}\n");

    base = {0 : "A", 1: "G", 2: "C", 3: "T"}
    rstr = ""
    for d in counts:
        if (max(d)/sum(d)) > cutoff:
            rstr += base[d.index(max(d))]
            # print(f"{base[d.index(max(d))]}", end="")
        else:
            rstr += "-"
            # print("-", end="")

    return rstr


def isfastq(f, verbose=False):
    """
    Simple check for fastq file
    """
    if 'fq' in f or 'fastq' in f:
        return True
    return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='input fastq file')
    parser.add_argument('-d', help='directory of fastq files')
    parser.add_argument('-k', help='kmer length (default=20)', type=int, default=20)
    parser.add_argument('-c', help="cutoff (default=0.5)", type=float, default=0.4)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    files = []
    if args.f:
        if os.path.exists(args.f):
            files.append(args.f)
    if args.d:
        files = files + [args.d + "/" + x for x in list(filter(isfastq, os.listdir(args.d)))]
    
    for f in files:
        c = count_kmers(f, args.k, args.v)
        s = predict(c, args.c, args.v)
        print(f"{f}\t{s}")

