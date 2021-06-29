"""
Merge the counts and abstracts. We take the abstracts file from partie, and the counts generated
by SearchSRA and create a single file with both all counts and average counts
"""

import os
import sys
import argparse
import gzip


__author__ = 'Rob Edwards'


class colors(object):
    

    color = {
        'HEADER': '\033[95m',
        'OKBLUE': '\033[94m',
        'OKGREEN': '\033[92m',
        'WARNING': '\033[93m',
        'FAIL': '\033[91m',
        'ENDC': '\033[0m',
        'BOLD': '\033[1m',
        'UNDERLINE': '\033[4m',
        'PINK': '\033[95m',
        'BLUE': '\033[94m',
        'GREEN': '\033[92m',
        'YELLOW': '\033[93m',
        'RED': '\033[91m',
        'WHITE': '\033[0m',
        }


def message(msg, color):
    """
    Print a message to stderr using color
    :param msg: the message to print
    :param color: the color to use
    :return: nothing
    """

    color = color.upper()
    if color not in colors.color:
        raise ColorNotFoundError(f"There is no color {color}")

    if os.fstat(0) == os.fstat(1):
        #  stderr is not redirected
        sys.stderr.write(f"{colors.color[color]}{msg}{colors.color['ENDC']}\n")
    else:
        sys.stderr.write(f"{msg}\n")

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
    if verbose:
        message(f"Reading counts_dir", "GREEN")
    for f in os.listdir(counts_dir):
        with open(os.path.join(counts_dir, f), 'r') as f:
            for l in f:
                if l.startswith('Sample'):
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
        message("Received {len(counts)} counts and {len(all_contigs)} contigs", "PINK")

    if verbose:
        message("Reading the abstracts and writing the output", "GREEN")

    if isinstance(all_contigs, set):
        all_contigs = sorted(all_contigs)

    if abstractsf.endswith('.gz'):
        abst = gzip.open(abstractsf, 'rt')
    else:
        abst = open(abstractsf, 'r')

    with open(reads_per_sample_file, 'w') as reads_out, open(average_per_proj, 'w') as average_out:
        average_out.write("Project\tTitle\tAbstract\tAnnotation\tComment\tTotal hits")
        for c in all_contigs:
            average_out.write(f"\t{c}")
        average_out.write("\n")

        for l in abst:
            # SRA Project     Title   Abstract        Annotation      Comment Runs
            p = l.strip().split("\t")
            if len(p) == 4 :
                project, title, abstract, runs = p
                annotation = ""
                comment = ""
            elif len(p) == 6:
                project, title, abstract, annotation, comment, runs = p
            else:
                if verbose:
                    message(f"Malformed Abstracts (len p: {len(p)}: {l}", "RED")
                continue



            runids = runs.split(',')
            run_counts = {c:[] for c in all_contigs}
            rowcount = 0
            for r in runids:
                if r not in counts:
                    message(f"Run {r} not found", "RED")
                    continue
                for c in all_contigs:
                    if c in counts[r]:
                        run_counts[c].append(counts[r][c])
                        rowcount += counts[r][c]
                        reads_out.write(f"{project}\t{r}\t{c}\t{counts[r][c]}\n")
            if rowcount > 0:
                average_out.write("\t".join([project, title, abstract, annotation, comment]))
                num = len(runids)
                average_out.write(f"\t{rowcount/num}")
                for c in all_contigs:
                    average_out.write(f"\t{sum(run_counts[c])/num}")
                average_out.write("\n")
    abst.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-d', help='directory of SRA Run counts per sample', nargs='+')
    parser.add_argument('-a', help='abstracts file that includes annotations and comments', required=True)
    parser.add_argument('-r', help='file to write reads per sample to', required=True)
    parser.add_argument('-p', help='file to write average counts per project to', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    countsperdir = {}
    contigsperdir = set()
    for subdir in args.d:
        cpr, ac = counts_per_sample(subdir, args.v)
        for c in cpr:
            if c not in countsperdir:
                countsperdir[c] = {}
            countsperdir[c].update(cpr[c])
        contigsperdir.update(ac)
    read_abstracts(args.a, args.r, args.p, countsperdir, contigsperdir, args.v)
