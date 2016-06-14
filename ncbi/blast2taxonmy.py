"""
Read a blast file that has gi| numbers and add the taxonomy to the end of the reads. We will add

kingdom phylum class order family genus species
"""

import os
import sys

import argparse
import taxon
import re

parser = argparse.ArgumentParser(description="Add taxonomy to a blast output file")
parser.add_argument('-b', help='blast ouput file', required=True)
parser.add_argument('-p', help='blast program (blastn, blastp, or blastx (default)', default='blastx')
args = parser.parse_args()

dbtype='nucl'
if args.p == 'blastx' or args.p == 'blastp':
    dbtype='prot'

sys.stderr.write("Reading taxonomy\n")
taxa=taxon.readNodes()
names,blastname = taxon.readNames()
divs = taxon.readDivisions()
gi2tax = taxon.readGiTaxId(dtype=dbtype)


want = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

with open(args.b, 'r') as f:
    for l in f:
        p=l.strip().split("\t")
        m = []
        if 'gi|' in p[0]:
            m = re.findall('gi\|(\d+)', p[0])
        elif 'gi|' in p[1]:
            m = re.findall('gi\|(\d+)', p[0])
        if m == []:
            continue
        if m[0] not in gi2tax:
            continue

        tid = gi2tax[m[0]]
        level = {}
        while taxa[tid].parent != '1' and tid != '1':
            if taxa[tid].rank in want:
                level[taxa[tid].rank] = names[tid].name
            tid = taxa[tid].parent

        for w in want:
            if w in level:
                p.append(level[w])
            else:
                p.append("")

        print("\t".join(p))




