"""
Calculate the mutation frequencies at each position from a multi fasta alignment
"""

import os
import sys

import argparse
from roblib import sequences

count = []
total = []
for i in range(100000):
    count.append({"A": 0, "G": 0, "T": 0, "C": 0, "N": 0})
    total.append(0)
longestseq = 0
bases = {"A", "G", "C", "T", "N"}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate the mutation frequencies at each position from a multi fasta alignment")
    parser.add_argument('-f', help='fasta alignment file')
    args = parser.parse_args()



    for (seqid, seq) in sequences.stream_fasta(args.f):
        if len(seq)  > longestseq:
            longestseq = len(seq)

        for i in range(len(seq)):
            if seq[i] in bases:
                count[i][seq[i]] += 1
                total[i] += 1

    print("A\tG\tC\tT\tN\tnumber of reads")
    for i in range(longestseq-1):
        if total[i] == 0:
            print("0\t0\t0\t0\t0\t0")
            continue
        results = []
        for b in "AGCTN":
            # results.append(round((1.0 * count[i][b] / total[i]), 2))
            results.append((1.0 * count[i][b] / total[i]))
        results.append(total[i])
        print("\t".join(map(str, results)))