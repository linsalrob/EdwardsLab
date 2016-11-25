"""
Create a substitution matrix based on Jims Rules

Current scoring approach

For every pair of bases in the alignment:
    if the base is the same score 0
    if the substitution is purine -> purine / pyrimidine -> pyrimidine (i.e. A <-> G or C <-> T) score 0.5
    if the substitution is purine -> pyrimidine (ie. {A, G} <-> {C, T}) score 1

Here we just calculate the score as a dict and optionally print it out

"""

import os
import sys
import argparse
from itertools import product

__author__ = 'Rob Edwards'


def score(length, alphabet={'A', 'T', 'G', 'C', '-'}):
    """
    Calculate the scoring matrix
    :param length: the length of the sequence to generate
    :param alphabet: the set of characters to use in the alphabet (default = {'A', 'T', 'G', 'C', '-'})
    :return: The dictionary of scores
    """

    scores = {}

    poss = list(product(alphabet, repeat=length)) # we need to convert this to a list so we can iterate it multiple times

    for i, tple in enumerate(poss):
        seq1 = "".join(tple)
        scores[seq1] = {}
        scores[seq1][seq1] = 0
        for j in range(i+1, len(poss)):
            seq2 = "".join(poss[j])
            scores[seq1][seq2] = 0
            for position in range(length):
                base1 = tple[position]
                base2 = poss[j][position]

                if base1 == base2:
                    # increment by 0
                    continue

                if (base1 == 'A' and base2 == 'T') or (base1 == 'T' and base2 == 'A'):
                    scores[seq1][seq2]+=0.5
                    continue

                if (base1 == 'G' and base2 == 'C') or (base1 == 'C' and base2 == 'G'):
                    scores[seq1][seq2] += 0.5
                    continue

                scores[seq1][seq2] += 1

    # now we just need to create the other half of the matrix
    for s in scores:
        for t in scores:
            if t not in scores[s]:
                scores[s][t] = scores[t][s]

    return scores


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a scoring matrix')
    parser.add_argument('-l', help='length of sequence to make scoring matrix for', required=True, type=int)
    args = parser.parse_args()

    sm = score(args.l)

    letters = list(sm.keys())
    letters.sort()

    print("\t" + "\t".join(letters))
    for l in letters:
        sys.stdout.write(l)
        for j in letters:
            sys.stdout.write("\t{}".format(sm[l][j]))
        sys.stdout.write("\n")