"""
Merge the counts and abstracts. We take the abstracts file from partie, and the counts generated
by SearchSRA and create a single file with both all counts and average counts
"""

import os
import sys
import argparse
import gzip

from .colours import message

def counts_per_sample(counts_dir, verbose=False):
    """
    Read the counts file and return a dict of the number of counts per sample
    We expect the counts file to have three columns, [sra_id, contig_id, count]. This is the output from
    process.smk using idxstats.

    :param str counts_dir: the directory of counts files
    :param bool verbose: more outputs
    :return dict[str, int]: the counts per sample, and a set of the contigs
    :rtype: dict[str, dict[str, int]], set
    """

    counts = {}
    all_contigs = set()
    for f in os.listdir(counts_dir):
        if verbose:
            message(f"Reading {os.path.join(counts_dir, f)}", "GREEN")
        with open(os.path.join(counts_dir, f), 'r') as f:
            for l in f:
                if l.startswith('sample'):
                    continue
                p = l.strip().split("\t")
                if int(p[2]) == 0:
                    continue
                if p[0] not in counts:
                    counts[p[0]] = {}
                counts[p[0]][p[1]] = int(p[2])
                all_contigs.add(p[1])

    return counts, all_contigs

def read_abstracts(abstractsf, reads_per_sample_file, average_per_proj, counts, all_contigs, verbose=False):
    """
    Read the abstracts file, and write the reads per sample and the average per project
    """

    if verbose:
        message("Reading the abstracts and writing the output", "GREEN")

    if isinstance(all_contigs, set):
        all_contigs = sorted(all_contigs)

    if abstractsf.endswith('.gz'):
        abst = gzip.open(abstractsf, 'rt')
    else:
        abst = open(abstractsf, 'r')

    with open(reads_per_sample_file, 'w') as reads_out, open(average_per_proj, 'w') as average_out:
        average_per_proj.write("Project\tTitle\tAbstract\tAnnotation\tComment")
        for c in all_contigs:
            average_per_proj.write(f"\t{c}")
        average_per_proj.write("\n")

        for l in abst:
            # SRA Project     Title   Abstract        Annotation      Comment Runs
            project, title, abstract, annotation, comment, runs = l.strip().split("\t")
            runids = runs.split(',')
            run_counts = {c:[] for c in all_contigs}
            for r in runids:
                if r not in counts:
                    message(f"Run {r} not found", "RED")
                    continue
                for c in all_contigs:
                    if c in counts[r]:
                        run_counts[c].append(counts[r][c])
                        reads_out.write(f"{project}\t{r}\t{c}\t{counts[r]}\n")
            average_per_proj.write("\t".join([project, title, abstract, annotation, comment]))
            num = len(runids)
            for c in all_contigs:
                average_per_proj.write(f"\t{sum(run_counts[c])/num}")
            average_per_proj.write("\n")
    abst.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-d', help='directory of SRA Run counts per sample', required=True)
    parser.add_argument('-a', help='abstracts file that includes annotations and comments', required=True)
    parser.add_argument('-r', help='file to write reads per sample to', required=True)
    parser.add_argument('-p', help='file to write average counts per project to', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    cpr, ac = counts_per_sample(args.d, args.v)
    read_abstracts(args.a, args.r, args.p, cpr, ac, args.v)