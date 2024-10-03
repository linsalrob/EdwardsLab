"""
Note: written by chatGPT  4o

This takes a piece of DNA and randomises the sequence while keeping the translation in frame one the same
"""

import os
import sys
import argparse
import random
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
from Bio.Data import CodonTable
__author__ = 'ChatGPT 4o'

def randomize_dna_sequence(dna_sequence):
    # Create a BioPython Seq object
    dna_seq = Seq(dna_sequence)
    # Get the translation table
    standard_table = CodonTable.unambiguous_dna_by_id[1]

    # Translate the original sequence
    protein_sequence = dna_seq.translate(to_stop=True)

    # Get all codons for each amino acid
    codon_dict = {aa: [] for aa in standard_table.protein_alphabet}
    for codon, aa in standard_table.forward_table.items():
        codon_dict[aa].append(codon)

    # Randomize the DNA sequence
    randomized_dna_sequence = ''
    for aa in protein_sequence:
        if aa == '*':  # Stop codon
            break
        codon = random.choice(codon_dict[aa])
        randomized_dna_sequence += codon

    return randomized_dna_sequence

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Randomise a DNA sequence, keeping the translationt the same')
    parser.add_argument('-d', '--dna', help='dna sequence')
    parser.add_argument('-v', '-verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not args.dna:
        parser.print_help()
        print("Try the sequence: ATGCGTATGCGTAA", file=sys.stderr)
        sys.exit()


    randomized_sequence = randomize_dna_sequence(args.dna)
    print(f"Original DNA sequence: {args.dna}")
    trans = Seq(args.dna).translate()
    print(f"Original translation: {trans}")

    print(f"Randomized DNA sequence: {randomized_sequence}")
    trans = Seq(randomized_sequence).translate()
    print(f"Original translation: {trans}")
