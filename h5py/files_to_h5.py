"""
Convert a series of files to h5 format.

Here we have a series of files

f1.tsv, f2.tsv, f3.tsv

where each has

[contig1, value]
[contig2, value]
[contig3, value]

and we want to create h5 objects.

There are two ways we could do this. For what I want (correlations) we should use contigs as the
groups in h5 and filenames to order the values.

Alternatively we could use filenames as groups and contigs to order the values
"""

import os
import sys
import argparse
import h5py
import numpy as np

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert a series of data to h5py format')
    parser.add_argument('-f', '--files', help='input files', nargs='+')
    parser.add_argument('-o', '--output', help='output h5 file to write', required=True)
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    parser.add_argument('-s', '--skip', help='skip lines that start *. This is very specific to parsed bam files',
                        action='store_true')
    args = parser.parse_args()

    data = {}
    all_contigs = set()
    for fname in args.files:
        with open(fname, 'r') as f:
            data[fname] = {}
            for l in f:
                if args.skip and l.startswith('*'):
                    continue
                p = l.strip().split("\t")
                data[fname][p[0]] = int(p[1])
                all_contigs.add(p[0])

    sorted_contigs = sorted(all_contigs)
    sorted_files = sorted(args.files)
    with h5py.File(args.output, "w") as f:
        grp = f.create_group("contigs")
        for contig in sorted_contigs:
            d = [data[f][contig] if contig in data[f] else 0 for f in sorted_files]
            dset = grp.create_dataset(contig, data=np.array(d))
