# -*- coding: UTF-8 -*-

"""
Score an alignment. Here is the general idea:

As part of moving the survival grant proposal in the limits direction, it seems like it would be great to create a
good analytical and visual metric for assessing the state of selection. Here is what seems like would be useful for
the proposal. If it has already been done let me know and I will come up with something else, but it isn’t easy.
Okay, take a sequence, e.g. 16S (it expands to meta-g later) compare it base by base to another closely related
sequence, score +1 for a single base like same substitution and -1 for switch substitution, i.e. purine to
pyrimidine, repeat along the whole sequence and get a score. To avoid averaging to zero all single substitutions
could be counted as +1. Repeat this with 2 base substitutions and again with 3 base substitutions. This gives 3 scores
for this first comparison. Repeat this for all 16S sequences. Take the triplet scores for each and plot them as
x, y and z on a graph. This defines a cloud of 16S scores that shows the density and distribution of the sequence
diversity. Realistically, the x values might combine the singles and doubles, the y values might combine the triples
and quadruples and the z values the quintuples and sextuples. In any case the idea is that a 16S cloud is defined.
From that cloud extract 3 values that describe deviation of the cloud from a sphere, i.e. prolation, oblation and
roughness. Now repeat this whole process for all of the sequences in the metagenome that can be linked to the 16S.
Get extracted deviation values for each gene set. These values are then used to plot a cloud of the entire gene set for
that individual (ideally), or more likely for that population. Extract the deviation values for that cloud to get
a description of the population shape. Then repeat this for each OUT to get a community cloud.
"""


import os
import sys
import argparse
from roblib import read_fasta

import substitution_rules

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='fasta file of sequence alignments', required=True)
    parser.add_argument('-s', help='maximum size of kmers to count (default=3)', default=3, type=int)
    args = parser.parse_args()

    fa = read_fasta(args.f)
    for sid in fa:
        fa[sid] = fa[sid].upper()

    sm = [""] # there is no point counting for size = 0, but it makes it easier if we are 0 indexed :)
    for i in range(1, args.s+1):
        sm.append(substitution_rules.score(i))

    # score all pairwise comparisons
    seqs = list(fa.keys())
    for i, seqid1 in enumerate(seqs):
        for j in range(i+1, len(seqs)):
            seqid2 = seqs[j]
            s1 = fa[seqid1]
            s2 = fa[seqid2]
            scores = []
            for sublen in range(1, args.s+1):
                score = 0
                posn = 0
                while posn + sublen < len(s1):
                    score += sm[sublen][s1[posn:posn+sublen]][s2[posn:posn+sublen]]
                    posn += 1
                scores.append(score)
            sys.stdout.write("{}\t{}\t".format(seqid1, seqid2))
            sys.stdout.write("\t".join(map(str, scores)))
            sys.stdout.write("\n")
            #print("{}\t{}\t{}\t{}\t{}".format(seqid1, seqid2, *scores))