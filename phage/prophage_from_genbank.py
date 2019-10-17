"""
Extract prophage regions from a genbank file
"""

import os
import sys
import argparse

from roblib import bcolors
from roblib import seqio_filter
from Bio import SeqIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-f', help='Genbank file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    input = seqio_filter.SeqioFilter(SeqIO.parse(args.f, "genbank"))

    amphage = False
    phagestart = 0
    phagecount = 0
    cdscount = 0
    for entry in input:
        for cds in entry.get_features('CDS'):
                is_phage = int(cds.qualifiers['is_phage'][0]) if 'is_phage' in cds.qualifiers else 0
                if is_phage:
                    cdscount += 1
                if is_phage and not amphage:
                    phagecount += 1
                    phagestart = cds.start
                    if cds.stop < phagestart:
                        phagestart = cds.stop
                    amphage = True
                if not is_phage and amphage:
                    amphage = False
                    # put the end of the phage at the start of this gene
                    phagestop = cds.stop
                    if cds.start < phagestop:
                        phagestop = cds.start
                    print(f"{args.f} | {phagecount} | {phagestart} | {phagestop} | {cdscount}")
                    cdscount = 0