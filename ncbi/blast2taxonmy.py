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
taxa=taxon.read_nodes()
names,blastname = taxon.read_names()
divs = taxon.read_divisions()
gi2tax = taxon.read_gi_tax_id(dtype=dbtype)

sys.stderr.write("Read taxonomy\n")

want = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

results = {}
with open(args.b, 'r') as f:
    for l in f:
        p=l.strip().split("\t")
        m = []
        if 'gi|' in p[0]:
            m = re.findall('gi\|(\d+)', p[0])
        elif 'gi|' in p[1]:
            m = re.findall('gi\|(\d+)', p[1])
        if m == []:
            continue

        if m[0] in results:
            p.append(results[m[0]])
            print("\t".join(p))
            continue

        if m[0] not in gi2tax:
            continue

        tid = gi2tax[m[0]]
        level = {}
        while tid != '0' and tid != '1' and tid in taxa and taxa[tid].parent != '1':
            if taxa[tid].rank in want:
                level[taxa[tid].rank] = names[tid].name
            tid = taxa[tid].parent

        resultstr = ""
        for w in want:
            if w != 'superkingdom':
                resultstr  += "\t"
            if w in level:
                p.append(level[w])
                resultstr  += level[w]
            else:
                p.append("")
                resultstr  += ""

        print("\t".join(p))




