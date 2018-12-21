"""
Calculate the k-mer bias for Kyle because he doesn't understand perl
"""

import os
import sys
import argparse
from itertools import combinations

def repN(n, original, sample, khash, verbose=False):
    """
    Permute the string, replacing each individual occurrence of N with {A, G, C, T}
    :param n: The test string
    :param original: The original string
    :param sample: The specific sample
    :param khash: our dict of dicts
    :param verbose: more output, but it is written to stderr
    :return:
    """

    if "N" not in n:
        return khash

    for base in ["A", "C", "G", "T"]:
        a = n.replace('N', base, 1)
        if "N" in a:
            repN(a, original, sample, khash, verbose)
        else:
            khash[sample][original] += khash[sample][a]
            if verbose:
                sys.stderr.write(f"{a}\t{original}\n")
    return khash


def read_kmer_counts(infile, verbose=False):
    """
    Read the k-mer counts file and create a dict of dicts
    :param infile: the file to read. This should be the output from GetNucFrequency_PerSeq_varK.pl
    :param verbose: more output, but it is written to stderr
    :return: a dict of dicts
    """

    khash = {}
    with open(infile, 'r') as f:
        names = f.readline().strip().split("\t")
        # define an empty set of dicts
        khash = {x:{} for x in names[1:]}
        for l in f:
            # convert this to an array
            p = l.strip().split("\t")
            key = p[0]
            for i in range(1, len(p)):
                khash[names[i]][key] = float(p[i])
    return khash


def process_all_kmers(kmers, khash, verbose=False):
    """
    Process and normalize all the k-mers
    :param kmers: an array of all the k-mers
    :param khash: our dict of dicts we've parsed out
    :param verbose: more output, but it is written to stderr
    :return:
    """

    samples = sorted(khash.keys())
    print("Matrix\t{}".format("\t".join(samples)))


    for k in kmers:
        sys.stdout.write("{}".format(k))
        for s in samples:
            if k not in khash[s]:
                sys.stdout.write("\t0")
                continue
            ksize = len(k)
            score = 1
            current = 0
            for i in range(ksize, 0, -1):
                for c in combinations(range(ksize), i):
                    testk = ""
                    startEnd = c[-1] - c[0] + 1
                    size = len(c)
                    if verbose:
                        sys.stderr.write(f"c: {c} Startend es {startEnd} y size {size}\n")
                    if startEnd > size:
                        # For each possible K-mer you must add the probabilities and that sum is the cell by which you will multiply.
                        count = c[0]
                        t = 0
                        numN = 0
                        while t < len(c):
                            if c[t] == count:
                                testk += k[c[t]] # append the letter
                                count += 1
                                t += 1
                            else:
                                testk += 'N'
                                count += 1
                                numN += 1
                        # Now you have to define the Freq value with N based on the sum of all the N components
                        if testk not in khash[s]:
                            khash[s][testk] = 0
                            khash = repN(testk, testk, s, khash, verbose)
                    else:
                        if verbose:
                            sys.stderr.write(f"Getting testk from: {c[0]} to {c[-1]} ")
                        for t in range(c[0], c[-1]+1):
                            testk += k[t]
                            if verbose:
                                sys.stderr.write(f"{t} ({testk}) ")
                        if verbose:
                            sys.stderr.write("\n")
                    if verbose:
                        sys.stderr.write(f"At this point: s: {s}  testk: {testk} score: {score} khash: {khash[s][testk]} current: {current}\n")
                    if 0 == current:
                        score *= khash[s][testk]
                    else:
                        score = score / khash[s][testk]
                    if verbose:
                        sys.stderr.write(f"At this point: s: {s}  testk: {testk} score: {score} khash: {khash[s][testk]} current: {current}\n")
                if 0 == current:
                    current = 1
                else:
                    current = 0
            sys.stdout.write(f"\t{score}")
        sys.stdout.write("\n")









__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Save Kyles butt')
    parser.add_argument('-f', help='output file from GetNucFrequency_PerSeq_varK.pl', required=True)
    parser.add_argument('-k', help='kmers to test. You can either specify multiple -k or use a comma separated list', required=True, action='append')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    khash = read_kmer_counts(args.f)

    # figure out all the kmers
    allkmers = set()
    for k in args.k:
        allkmers.update(k.split(","))

    process_all_kmers(allkmers, khash, args.v)
