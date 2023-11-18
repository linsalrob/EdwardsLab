"""
We use the filter metagenomes script but do it recursively from 1 to 2000 reads matching
"""
import gzip
import os
import sys
import argparse
from pathlib import Path

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
    for path in Path(args.directory).rglob('*.tophit_aln.gz'):
        metagenome = os.path.split(path)[-2]
        count[metagenome] = {}
        with gzip.open(path, 'rt') as f:
            for l in f:
                p = l.strip().split("\t")
                if float(p[2]) < args.percent:
                    continue
                count[metagenome][p[0]] = count[metagenome].get(p[0], 0) + 1

    print("Min number of proteins\tNumber of metagenomes")
    for i in range(2000):
        mgxcount = 0
        for m in count:
            if len(count[m].keys()) > i:
                mgxcount += 1
        print(f"{i}\t{mgxcount}")
