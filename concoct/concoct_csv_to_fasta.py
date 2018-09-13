"""
Convert the output from concoct csv to a set of fasta files one per bin
"""

import os
import sys
import argparse
from roblib import stream_fasta

__author__ = 'Rob Edwards'

def read_csv(csvf, verbose=False):
    """
    Read the csv list and return a hash of contigs->bins and the number of bins
    :param csvf: csv file to read
    :param verbose: more output
    :return: a dict of contigs -> bins and an int of the largest bin number
    """
    bins={}
    x=0
    with open(csvf, 'r') as f:
        for l in f:
            p = l.strip().split(",")
            bins[p[0]]=int(p[1])
            if int(p[1]) > x:
                x = int(p[1])
    return bins, x


def write_fasta_files(faf, odir, bins, maxb, verbose=False):
    """
    Read the sequences from faf and write them into a set of files in odir.
    :param faf: The source fasta file
    :param odir: the output directory
    :param bins: the hash of contigs -> bin
    :param maxb: the maximum bin number
    :param verbose: more output
    :return: nada
    """

    if not os.path.exists(odir):
        os.mkdir(odir)
    
    outputfiles = []
    for i in range(maxb+1):
        outputfiles.append(open(os.path.join(odir, f"bin_{i}.fna"), 'w'))

    written_to=set()

    for fa, seq in stream_fasta(faf, True):
        faid = fa.split(" ")[0]
        if faid not in bins:
            if verbose:
                sys.stderr.write(f"Sequence {faid} not found in a bin\n")
            continue
        outputfiles[bins[faid]].write(">{}\n{}\n".format(fa, seq))
        written_to.add(bins[faid])

    for o in outputfiles:
        o.close()

    for i in range(maxb+1):
        if i not in written_to:
            os.remove(os.path.join(odir, f"bin_{i}.fna"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert the csv contig id,bin to a set of fasta files')
    parser.add_argument('-f', help='fasta file of DNA sequences', required=True)
    parser.add_argument('-c', help='concoct csv file', required=True)
    parser.add_argument('-d', help='output directory', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    bins, maxb = read_csv(args.c, args.v)
    write_fasta_files(args.f, args.d, bins, maxb, args.v)
