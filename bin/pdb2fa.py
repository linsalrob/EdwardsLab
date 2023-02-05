"""
Convert a PDB file to a fasta file

Taken from https://www.biostars.org/p/435629/
"""

import os
import sys
from Bio import SeqIO
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert a PDB file to a fasta file')
    parser.add_argument('-p', help='input PDB file', required=True)
    parser.add_argument('-f', help='output fasta file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    with open(args.p, 'r') as pdb_file, open(args.f, 'w') as fasta_file:
        for record in SeqIO.parse(pdb_file, 'pdb-atom'):
            if record.id.startswith("???"):
                print(f">{args.p}\n{record.seq}", file=fasta_file)
            else:
                print(f">{record.id}\n{record.seq}", file=fasta_file)


