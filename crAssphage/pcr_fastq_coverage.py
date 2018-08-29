"""
Given a directory of files, figure out the average coverage of our PCR products.
"""

import os
import sys
import argparse
import re
from roblib import stream_fastq

def parse_dir(dir, verbose=False):
    """
    Parse the directory of files
    :param dir:
    :param verbose:
    :return:
    """

    lengths = {}
    for f in os.listdir(dir):
        m = re.search('^(\w+)_(Primer\w)_', f)
        if not m:
            sys.stderr.write("Error: can't parse {}\n".format(f))
            continue
        (srr, primer) = m.groups()
        if srr not in lengths:
            lengths[srr] = {'PrimerA' : 0, 'PrimerB' : 0, 'PrimerC' : 0}
        for seqid, header, seq, qualscores in stream_fastq(os.path.join(dir, f)):
            lengths[srr][primer] += len(seq)
    return lengths

def print_coverage(lengths):
    """
    Print the coverage of the amplicons
    :param lengths: The dict of SRR and lengths
    :return:
    """

    primerAlen = 1330
    primerBlen = 1353
    primerClen = 1237

    print("SRR ID\tPrimer A\tPrimer B\tPrimer C")
    for s in lengths:
        sys.stdout.write(s)
        sys.stdout.write("\t{}".format(1.0 * lengths[s]["PrimerA"]/primerAlen))
        sys.stdout.write("\t{}".format(1.0 * lengths[s]["PrimerB"]/primerBlen))
        sys.stdout.write("\t{}\n".format(1.0 * lengths[s]["PrimerC"]/primerClen))





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-d', help='directory of files', required=True)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    lens = parse_dir(args.d, args.v)
    print_coverage(lens)