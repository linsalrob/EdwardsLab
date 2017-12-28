"""
Read the best hit, get the taxa of that hit and report the number of [kingdom/phyla/genera/species] vs the number of proteins
"""

import os
import sys
import argparse
import taxon
import re

ffdir = 'flatfiles' # the location of our genbank flatfiles
blastpf = 'phage_genes.nr.blastp'

genome = {} # hash with protein IDs as key and genome as value
numprots = {} # hash with genome ID as key and number of proteins as value

taxondir = "/home2/db/taxonomy/current/"
sys.stderr.write("Reading taxonomy\n")
taxa = taxon.read_nodes(directory=taxondir)
names, blastname = taxon.read_names(directory=taxondir)
sys.stderr.write("Read taxonomy\n")

# first read how many proteins there are per genome
# and create a list of protein ids->genomes.
# Note also throw a fatal error if IDs duplicated because I forgot to check this earlier!
for f in os.listdir(ffdir):
    with open(os.path.join(ffdir, f), 'r') as fin:
        for l in fin:
            p=l.strip().split("\t")
            if p[5] in genome:
                sys.stderr.write("FATAL AND BUGGER: {} was found in {} and {}\n".format(p[5], genome[p[0]], p[0]))
            genome[p[5]] = p[0]
            numprots[p[0]] = numprots.get(p[0], 0) + 1

prottaxa = {} # hash of protein taxa
want = ['superkingdom', 'phylum', 'genus', 'species']
with open(blastpf, 'r') as fin:
    for l in fin:
        p = l.strip().split("\t")
        if p[0] in prottaxa:
            continue
        for tid in p[14].split(";"):
            level = {}

            while tid != '0' and tid != '1' and tid in taxa and taxa[tid].parent != '1':
                if taxa[tid].rank in want:
                    level[taxa[tid].rank] = names[tid].name
                tid = taxa[tid].parent
            # find the first entry that has all the taxonomic levels we want
            keep = True
            for w in want:
                if w not in level:
                    keep = False
            if keep:
                prottaxa[p[0]] = level
                break

hostlocation={}
with open('phage_host_location.txt', 'r') as fin:
    for l in fin:
        p=l.strip().split("\t")
        hostlocation[p[0]] = [p[2], p[3]]

hitcount = {}
genometaxa = {}
for p in prottaxa:
    if genome[p] not in genometaxa:
        genometaxa[genome[p]] = {}
        hitcount[genome[p]] = 0
        for w in want:
            genometaxa[genome[p]][w] = set()
    for w in prottaxa[p]:
        genometaxa[genome[p]][w].add(prottaxa[p][w])
    hitcount[genome[p]] += 1

print("#Genome ID\tNum proteins\tNum proteins with hits\tHost\tBody location\tKingdom\tPhylum\tGenus\tSpecies")

for g in genometaxa:
    sys.stdout.write("{}\t{}\t{}\t".format(g, numprots[g], hitcount[g]))
    if g in hostlocation:
        sys.stdout.write("\t".join(hostlocation[g]))
    else:
        sys.stdout.write("-\t-\t")
    for w in want:
        sys.stdout.write("\t{}".format(len(genometaxa[g][w])))
    sys.stdout.write("\n")


