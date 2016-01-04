"""
Generate a collectors curve of kmer profiles from a directory of alignment files

We assume that all alignments start at the same point. This is probably true but you can check it using
perl -ne 'if (substr($_, 0, 1) !~ /[>-]/) {print substr($_, 0, 10), "\n"}' all.aln | less

We also assume that alignments are all the same length, and you can check tht with
perl -ne 'if ($_ !~ /^>/) {print length($_), "\n"}' all.aln | sort -u


"""

import os
import sys
import argparse

from roblib import sequences


def average_genotype_frequency(seqs, all_genotypes, kmer=10, min_num_reads=1, verbose=False):
    """
    Calculate the average genotype frequency using k-mer profiles. At each position
    genotypes must be in >= min_num_reads to be included

    :param verbose: Print extra output
    :type verbose: bool
    :param all_genotypes: List of set of all genotypes at a position in the alignment. The list must be the same length
        as the first sequence in seqs
    :type all_genotypes: list of set
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

    if len(all_genotypes) != len(seqs[0]):
        sys.exit("The all_genotypes list (len {}) must be the same as the length of the sequences (len {})".format(
            len(all_genotypes), len(seqs[0])))

    sites = 0

    bases = {'A', 'T', 'G', 'C', 'a', 't', 'g', 'c', }

    for i in range(len(seqs[0]) - kmer):
        counts = {}
        for s in seqs:
            kseq = s[i:i + kmer]
            keep_kmer = True
            for k in kseq:
                if k not in bases:
                    keep_kmer = False
            if not keep_kmer:
                continue
            counts[kseq] = counts.get(kseq, 0) + 1

        for kseq in counts:
            if counts[kseq] >= min_num_reads:
                all_genotypes[i].add(kseq)
                sites += 1


    av = []
    for i in all_genotypes:
        sg = 0
        ng = 0
        if i:
            sg += len(i)
            ng += 1
        if ng > 0:
            av.append(1.0 * sg/ng)

    if verbose:
        sys.stderr.write("Found an additional {} sites and {} genotypes\n".format(sites, 1.0 * sum(av)/len(av)))

    return all_genotypes

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate a collectors curve of kmer profiles from a directory of alignment files")
    parser.add_argument('-d', help='directory of files to use')
    args = parser.parse_args()

    existing = []
    current_file_no = 0

    for aln_file in os.listdir(args.d):
        # sys.stderr.write("{}\t{}\n".format(current_file_no, aln_file))
        current_file_no += 1
        seqs = []
        seqids = []
        for sid, sq in sequences.stream_fasta(os.path.join(args.d, aln_file)):
            seqids.append(sid)
            seqs.append(sq)
        newresults = [set() for x in seqs[0]]
        if not existing:
            existing = [set() for x in seqs[0]]
        newresults = average_genotype_frequency(seqs, newresults, 10, 2, False)
        new_genotypes = []
        for i in range(len(newresults)):
            new = 0
            for kmer in newresults[i]:
                if kmer not in existing[i]:
                    existing[i].add(kmer)
                    new += 1
            if new > 0:
                new_genotypes.append(new)
        if len(new_genotypes) > 0:
            print("{}\t{}".format(current_file_no, 1.0 * sum(new_genotypes)/len(new_genotypes)))
        else:
            print("{}\t0".format(current_file_no))


