"""
Filter the sequences based on kmer and 16S composition
"""

import os
import sys

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot the SRA data from PARTIE")
    parser.add_argument('-p', help='partie output file', required=True)
    parser.add_argument('-e', help='Experiment library file (run id\texperiment library', required=True)
    parser.add_argument('-k', help='minimum value for k-mer composition to keep (0-1)', default=-10, type=float)
    parser.add_argument('-s', help='maximum value for 16S composition to keep (0-100)', default=200, type=float)
    parser.add_argument('-o', help='output a few of the ones that we want to keep that are amplicons', action='store_true')
    args = parser.parse_args()

    explib = {}
    with open(args.e, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            if len(p) == 1:
                p.append('OTHER')
            explib[p[0]] = p[1]

    data = {}
    experimentlibraries = {}
    with open(args.p, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            # data is unique kmers, percent 16S, percent phage, percent prok, percent prok + phage
            p[0] = p[0].replace('.sra', '')
            data[p[0]] = [float(p[2]), float(p[4]), float(p[6]), float(p[8]), float(p[6]) + float(p[8])]
            if p[0] in explib:
                if explib[p[0]] in experimentlibraries:
                    experimentlibraries[explib[p[0]]].append(p[0])
                else:
                    experimentlibraries[explib[p[0]]] = [p[0]]
            else:
                sys.stderr.write("No " + p[0] + " in exp\n")

    keep = set()

    printout=0
    for s in data.keys():
        if data[s][0] > args.k and data[s][1] < args.s:
            keep.add(s)
            if args.o and explib[s] == 'AMPLICON' and printout < 100:
                printout+=1
                print(s + "\tAMPLICON")


    # summarize how many we keep
    count = {}
    for s in keep:
        count[explib[s]] = count.get(explib[s], 0) + 1

    for c in count:
        print(c + "\t" + str(count[c]))
