"""
Extract genotypes from metagenomes by concatenating matching sequences.
"""

import sys

import argparse
from roblib import sequences
import re
import matplotlib.pyplot as plt


def replace_leading_trailing(seq):
    """
    Replace leading and trailing Ns with -s

    :param seq: the sequence
    :type seq: str
    :return: the sequence with leading trailing N's replaced
    :rtype: str
    """

    validbase = {"A", "G", "C", "T"}

    lastbase = 0
    inseq = False
    newseq = []
    for i in range(len(seq)):
        newseq.append("")
    for (i, j) in enumerate(seq):
        if j in validbase:
            inseq = True
            newseq[i] = j
            lastbase = i
        elif inseq:
            newseq[i] = j
        elif j == "N":
            newseq[i] = "-"
        else:
            newseq[i] = j

    newnewseq = newseq[0:lastbase]
    for i in range(lastbase+1, len(newseq)):
        if newseq[i] == "N":
            newnewseq.append("-")
        else:
            newnewseq.append(newseq[i])

    return "".join(newnewseq)


# data is going to be an array where each element is a hash
# the hash has to have the following elements:
#         seq :  this is the sequence and will be mutable - we will change the seqeunce
#         ids : this is a list of seqid that contribute to this sequence


count=[]
total = []
# we initially set this up to be 10kb. Hopefully no sequences longer than that!
for i in range(10000):
    count.append({"A": 0, "G": 0, "T": 0, "C": 0, "N": 0})
    total.append(0)
longestseq = 0
firstrun = True
bases = {'A', 'T', 'G', 'C', 'N'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extract genotypes from metagenomes by concatenating matching sequences")
    parser.add_argument('-f', help='fasta sequence alignment file', required=True)
    parser.add_argument('-n', help='minimum number of sequences a genotype must be in (default = 1)', default=1, type=int)
    parser.add_argument('-b', help='plot start position', type=int)
    parser.add_argument('-e', help='plot end position', type=int)
    parser.add_argument('-c', help='cutoff to print the values (e.g. 0.8)', type=float)
    parser.add_argument('-p', help='make the plot', action="store_true")
    parser.add_argument('-a', help='print all results by base (the default is to sort the possibilities)', action="store_true")

    args = parser.parse_args()


    # to start we are just going to merge identical sequences
    byseq = {}
    for (seqid, seq) in sequences.stream_fasta(args.f):
        seq = seq.upper() # make sure we have only upper case
        seq = seq.strip() # strip of leading/trailing whitespace
        #seq = seq.replace('N', '-')  # remove any N's in the sequence
        seq = replace_leading_trailing(seq)
        seq = seq.rstrip('-')  # remove all the trailing -s

        keep_seq = True
        for k in seq:
            if k not in bases and k != "-":
                if keep_seq:
                    sys.stderr.write("Skipped {} becuase it has base {}\n".format(seqid, k))
                keep_seq = False
        if not keep_seq:
            continue

        if len(seq) > longestseq:
            longestseq = len(seq)

        for (i, j) in enumerate(seq):
            if j in bases:
                count[i][j] = count[i][j] + 1
                total[i] += 1

    # now lets just make a plot of the data

    startseq = 0
    if args.e:
        longestseq = args.e
    if args.b:
        startseq = args.b

    if args.p:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    x = range(startseq, longestseq)

    toprint = []
    for base in "A", "T", "G", "C", "N":
        counts = []
        for i in range(startseq,longestseq):
            if not total[i]:
                continue
            if args.c and 1.0 * count[i][base] / total[i] < args.c and 1 - (1.0 * count[i][base] / total[i]) < args.c:
                toprint.append(i)
            counts.append(1.0 * count[i][base]/total[i])
        if args.p:
            ax.plot(x, counts, label=base)

    if args.c and toprint and not args.a:
        for p in toprint:
            results = []
            for base in "A", "T", "G", "C", "N":
                results.append(1.0 * count[p][base] / total[p])
            results.sort(reverse=True)
            print("\t".join(map(str, results)))


    if args.a and toprint:
        print("A\tG\tT\tC\tN\n")
        for p in toprint:
            print("{}\t{}\t{}\t{}\t{}\t{}".format(1.0 * count[p]["A"] / total[p],
                                          1.0 * count[p]["G"] / total[p],
                                          1.0 * count[p]["T"] / total[p],
                                          1.0 * count[p]["C"] / total[p],
                                          1.0 * count[p]["N"] / total[p],
                                              total[p]))

    if args.p:
        ax.legend()
        fig.set_facecolor('white')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()

        plt.show()

