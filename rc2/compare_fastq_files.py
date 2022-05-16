"""
Compare all the fastq files we have downloaded.

Note we use blake2b rather than md5sum because it is much faster:

See also: https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file
"""

import os
import sys
import argparse
import hashlib

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare all fastq files in a set of subdirectories')
    parser.add_argument('-d', help='super directory', required=True)
    parser.add_argument('-o', help='output file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    hashes = {}
    alldirs = []
    for f in os.listdir(args.d):
        subdir = os.path.join(args.d, f)
        if os.path.isdir(subdir):
            if args.v:
                sys.stderr.write(f"Reading {f}\n")
            alldirs.append(f)
            for fq in os.listdir(subdir):
                if fq.startswith('index'):
                    continue
                if not fq.endswith('.gz'):
                    sys.stderr.write(f"Not gzip compressed: {os.path.join(subdir, fq)}\n")
                    continue
                if fq not in hashes:
                    hashes[fq] = {}
                file_hash = hashlib.blake2b()
                with open(os.path.join(subdir, fq), 'rb') as fin:
                    while chunk := fin.read(8192):
                        file_hash.update(chunk)
                hashes[fq][f] = file_hash.hexdigest()

    if args.v:
        sys.stderr.write(f"Output to {args.o}\n")
    with open(args.o, 'w') as out:
        out.write("\t".join(["fastq"] + alldirs))
        out.write("\n")
        for fq in hashes:
            out.write(fq)
            same = set()
            for d in alldirs:
                if d in hashes[fq]:
                    out.write(f"\t{hashes[fq][d]}")
                    same.add(hashes[fq][d])
                else:
                    out.write("\t")
            if len(same) == 1:
                out.write("\tUnique\n")
            else:
                out.write("\tDifferent\n")
