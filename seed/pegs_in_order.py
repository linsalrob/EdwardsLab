"""
Write a list of the pegs in the order that they appear on the chromosome. We parse Features/peg/tbl for this information
"""

import os
import sys
import argparse
import re
import operator

parser = argparse.ArgumentParser(description="Generate a sorted list of pegs from Features/peg/tbl")
parser.add_argument('-d', help="SEED directory to generate pegs_in_order for", required=True)
args = parser.parse_args()

if not os.path.exists(os.path.join(args.d, 'Features/peg/tbl')):
    sys.exit("ERROR: {} not found\n".format(os.path.join(args.d, 'Features/peg/tbl')))

pegs = []

with open(os.path.join(args.d, 'Features/peg/tbl'), 'r') as f:
    for l in f:
        peg, loc=l.strip().split("\t")
        m = re.search('^(.*)_(\d+)_(\d+)$', loc)
        contig, start, stop = m.groups()
        pegs.append([peg, contig, start, stop])

pegs_sorted = sorted(pegs, key = operator.itemgetter(1,2,3))

for p in pegs_sorted:
    l = "_".join(map(str, p[1,2,3]))
    print("{}\t{}".format(p, l))

