"""
We use the filter metagenomes script but do it recursively from 1 to 2000 reads matching
"""
import gzip
import os
import sys
import argparse
from pathlib import Path
import re

def normalise():
    """
    Find out if we can normalise the data
    """

    normalized_data = {}
    if not os.path.exists("fastq_counts.tsv"):
        print("WARNING: We did not find a fastq_counts.tsv file so not normalizing to the number of reads", file=sys.stderr)
        return normalized_data 

    s = re.compile('/(\w+)_R1')
    with open("fastq_counts.tsv", "r") as f:
        for l in f:
            if "R1" in l:
                p = l.strip().split("\t")
                mg = s.findall(p[0])[0]
                normalized_data[mg] = int(p[1])
    return normalized_data

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-d', '--directory', help='base directory with the alignment files', required=True)
    parser.add_argument('-p', '--percent', help='minimum percent identity for a match (default: 90%)',
                        type=float, default=0.9)
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    if args.percent > 1:
        args.percent /= 100

    # find the files
    count = {}
    total = {}
    for path in Path(args.directory).rglob('*.tophit_aln.gz'):
        metagenome = os.path.split(path.parent)[1]
        count[metagenome] = {}
        total[metagenome] = 0
        with gzip.open(path, 'rt') as f:
            for l in f:
                p = l.strip().split("\t")
                if float(p[2]) < args.percent:
                    continue
                count[metagenome][p[0]] = count[metagenome].get(p[0], 0) + 1
                total[metagenome]+=1

    with open("cov_counts.tsv", 'w') as out:
        print("Min number of proteins\tNumber of metagenomes", file=out)
        for i in range(2000):
            mgxcount = 0
            for m in count:
                if len(count[m].keys()) > i:
                    mgxcount += 1
            print(f"{i}\t{mgxcount}", file=out)

    norm = normalise()
    with open("mg_read_frac.tsv", "w") as out:
        print("Metagenome\tNormalized number of reads", file=out)
        for m in total.keys():
            if m in norm:
                print(f"{m}\t{total[m]/norm[m]}", file=out)
            else:
                print(f"Error: no normalization: {m}", file=sys.stderr)
                print(f"{m}\t{total[m]}", file=out)

