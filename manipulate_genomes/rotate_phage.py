"""
Rotate a fasta sequence so that the start is just before the terminase
"""

import os
import sys
import argparse
from roblib import stream_blast_results, rc, read_fasta

def read_blast_file(bf, evalue, verbose=False):
    """
    Read the blast file and store tuples of contig, start, stop for each hit
    :param bf: blast file to parse
    :param evalue: E value to keep matches
    :param verbose: more output
    :return:
    """

    results = {}

    for br in stream_blast_results(bf, verbose=verbose):
        if br.evalue > evalue:
            continue
        if br.query not in results:
            results[br.query] = { '-ve' : set(), '+ve' : set() }
        strand = '+ve'
        (s, e) = (br.query_start, br.query_end)
        if br.query_end < br.query_start:
            strand = '-ve'
            (s, e) = (br.query_end, br.query_start)
            
        results[br.query][strand] = results[br.query][strand] | set(range(s, e))
    return results


def is_consecutive(l):
    """
    Determines whether a list contains consecutive numbers. If it does, the sum of the
    numbers should be the same as n * (min + max) / 2 (where n is num elements)
    :param l: A set of numbers
    :return: True or False
    """

    total = sum(l)
    mathsum = len(l) * (min(l) + max(l)) / 2
    return total == mathsum

def find_gene(brs):
    """
    Test for where the genes are
    :param brs: The blast results
    :return:
    """

    contigs = brs.keys()
    if len(contigs) == 1:
        sys.stderr.write("There was only one contig with a match ({})\n".format(list(contigs)[0]))
    else:
        sys.stderr.write("Multiple contigs had matches: {}\n".format(contigs))

    for c in contigs:
        sys.stderr.write(f"{c}\n")
        for strand in ['-ve', '+ve']:
            if len(brs[c][strand]) == 0:
                continue
            if is_consecutive(brs[c][strand]):
                sys.stderr.write("There is a consecutive match to {} {} {}\n".format(c, max(brs[c][strand]), min(brs[c][strand])))
            else:
                sys.stderr.write("There are multiple discontinuous matches to {}. Try adjusting the evalue\n".format(c))





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Recircularize the genome at the terminase")
    parser.add_argument('-f', help='fasta file of the genome', required=True)
    parser.add_argument('-l', help='blast output compared to terminase LARGE gene', required=True)
    parser.add_argument('-s', help='blast output compared to terminase SMALL gene', required=True)
    parser.add_argument('-e', help='evalue cutoff (default=all)', type=int, default=100)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    sys.stderr.write("Parsing small subunit\n")
    ss = read_blast_file(args.s, args.e, args.v)
    find_gene(ss)

    sys.stderr.write("Parsing large subunit\n")
    ls = read_blast_file(args.l, args.e, args.v)
    find_gene(ls)
