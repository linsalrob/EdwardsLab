"""
Extract the locus information, gene product, and translation from a genbank file
"""

import os
import sys

import argparse
from Bio import SeqIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extract the locus information, gene product, and translation from a genbank file")
    parser.add_argument('-f', help='genbank file', required=True)
    args = parser.parse_args()

    for seq in SeqIO.parse(args.f, 'genbank'):
        for feature in seq.features:
            pi = 'None'
            if 'protein_id' in feature.qualifiers:
                pi = feature.qualifiers['protein_id'][0]

            gs = "None"
            if 'gene' in feature.qualifiers:
                gs = feature.qualifiers['gene'][0]

            pd = 'None'
            if 'product' in feature.qualifiers:
                pd = feature.qualifiers['product'][0]

            tl = "None"
            if 'translation' in feature.qualifiers:
                tl = feature.qualifiers['translation'][0]

            if 'gpA' in gs or 'gpA' in pd:
                print("\t".join([pi, gs, pd, tl]))
