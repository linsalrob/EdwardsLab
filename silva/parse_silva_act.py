"""
This script parses the output from Silva ACT (https://www.arb-silva.de/aligner)

If you include Search and classify and then download the fasta file, you can parse that using this script
"""

import os
import sys
import argparse
import re
from roblib import stream_fasta

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse SILVA ACT')
    parser.add_argument('-f', help='input ACT fasta file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    s = re.compile(r'\[(.*?)\]')
    want = [
            'rdp_genus', 'lca_tax_rdp',
            'slv_genus', 'lca_tax_slv',
            'gtdb_genus', 'lca_tax_gtdb'
            ]

    for seqid, seq in stream_fasta(args.f, True):
        m = s.findall(seqid)
        if not m:
            print(f"Can't parse {seqid}", file=sys.stderr)
            continue
        result = ["", "", "", "", "", "", seqid, seq.replace('-', '')]
        for t in m:
            if '=' in t:
                p = t.split("=")
                who = p[0]
                what = p[1]
                last = ""
                if ';' in what:
                    last = what.split(';')[-2]
                if who in want:
                    idx = want.index(who)
                    result[idx-1] = last
                    result[idx] = what
        print("\t".join(result))
