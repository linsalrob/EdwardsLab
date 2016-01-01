"""
Calculate the average SNP frequency at all residues in the alignment. Assumes the alignment is in fasta format!
"""

import os
import sys
import argparse

from roblib import sequences


def average_snp_frequency(seqs, min_num_reads=1):
    """
    Calculate the average SNP frequency at sites where there is a SNP across the alignment. At each position
    SNPs must be in >= min_num_reads to be included

    :param sequences: An array of just the sequences. All elements in the array should be the same length (ie.
        sequences should be padded with -)
    :type sequences: list
    :param min_num_reads:minimum number of reads for a SNP to be in to be included
    :type min_num_reads: int
    :return:
    :rtype:
    """

    allsnps = []
    sites = 0

    bases = {'A', 'T', 'G', 'C', 'a', 't', 'g', 'c',}

    for i in range(len(seqs[0])):
        counts = {}
        for s in seqs:
            if s[i] not in bases:
                continue
            counts[s[i]] = counts.get(s[i], 0) + 1
        snps = 0
        for c in counts:
            if counts[c] >= min_num_reads:
                snps += 1
        if snps > 1:
            allsnps.append(snps)
            sites += 1

    sys.stdout.write("Alignment length: {} Informative sites: {} ".format(len(seqs[0]), sites))
    if sites == 0:
        # if there are no informative sites everything has the same sequence so we have one SNP!
        print("Average number of SNPs: 1")
    else:
       print("Average number of SNPs: {}".format(1.0 * sum(allsnps)/len(allsnps)))


def average_genotype_frequency(seqs, kmer=10, min_num_reads=1):
    """
    Calculate the average genotype frequency using k-mer profiles. At each position
    genotypes must be in >= min_num_reads to be included

    :param kmer: The length of the sequence to use to calculate the genotypes
    :type kmer: int
    :param sequences: An array of just the sequences. All elements in the array should be the same length (ie.
        sequences should be padded with -)
    :type sequences: list
    :param min_num_reads:minimum number of reads for a SNP to be in to be included
    :type min_num_reads: int
    :return:
    :rtype:
    """

    all_genotypes = []
    sites = 0

    bases = {'A', 'T', 'G', 'C', 'a', 't', 'g', 'c', }

    for i in range(len(seqs[0])-kmer):
        counts = {}
        for s in seqs:
            kseq = s[i:i+kmer]
            keep_kmer = True
            for k in kseq:
                if k not in bases:
                    keep_kmer = False
            if not keep_kmer:
                continue
            counts[kseq] = counts.get(kseq, 0) + 1
        snps = 0
        for c in counts:
            if counts[c] >= min_num_reads:
                snps += 1
        if snps > 1:
            all_genotypes.append(snps)
            sites += 1

    sys.stdout.write("Alignment length: {} Informative kmers: {} ".format(len(seqs[0]), sites))
    if sites == 0:
        print("Average number of genotypes: 1")
    else:
        print("Average number of genotypes: {}".format(1.0 * sum(all_genotypes) / len(all_genotypes)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate the average number of SNPs at informative sites")
    parser.add_argument('-a', help='Alignment file in fasta format')
    parser.add_argument('-s', help='Count SNPs', action='store_true')
    parser.add_argument('-g', help='Count genotypes', action='store_true')
    parser.add_argument('-k', help='k-mer size to estimate genotypes. Default=10', default=10, type=int)
    parser.add_argument('-m', help='minimum number of reads for site to be informative', default=2, type=int)
    args = parser.parse_args()

    seqids = []
    seqs = []

    for rid, seq in sequences.stream_fasta(args.a):
        seqids.append(rid)
        seqs.append(seq)

    if args.s:
        average_snp_frequency(seqs, min_num_reads=args.m)
    if args.g:
        average_genotype_frequency(seqs, kmer=args.k, min_num_reads=args.m)