"""
Test a genbank file and make sure all the protein_ids are unique
"""

import os
import sys
import argparse
from Bio import SeqIO

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-f', help='genbank file', required=True)
    args = parser.parse_args()

    pids = set()
    rc = 0
    for seq in SeqIO.parse(args.f, "genbank"):
        rc+=1;
        print(f"record {rc}: {seq.id}")
        for feat in seq.features:
            if feat.type != "CDS":
                continue
            if 'protein_id' not in feat.qualifiers:
                thisid = " ".join(feat.qualifiers.get('locus_tag', [str(feat.location)]))
                print(f"No protein id in {thisid}")
                continue
            pid = "|".join(feat.qualifiers["protein_id"])
            if pid in pids:
                print(f"{pid} is not unique")
            pids.add(pid)