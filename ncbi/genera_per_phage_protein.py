"""
Calculate the kingdom / phylum / genus / species per genera for the phages
"""

import os
import sys
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate the kingdom / phylum / genus / species per genera for the phages")
    parser.add_argument('-d', help='directory with phage flat files, one file per phage', required=True)
    parser.add_argument('-i', help='file with id, taxid, taxonomy (just kingdom / phylum / genus / species). Output from blast_tax_to_genera.py', required=True)
    parser.add_argument('-l', help='file with location in body (default: phage_host_location.txt)', default='phage_host_location.txt')
    parser.add_argument('-b', help='Only print phages for which we have a body site associated with the host', action='store_true')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    bodysite={}
    with open(args.l, 'r') as fin:
        for l in fin:
            p=l.strip().split("\t")
            bodysite[p[0]] = p[3]


    genome = {} # this is a hash of proteins -> genomes
    count = {}
    for f in os.listdir(args.d):
        if args.v:
            sys.stderr.write("Reading genome {}\n".format(f))
        with open(os.path.join(args.d, f), 'r') as fin:
            for l in fin:
                p=l.strip().split("\t")
                genome[p[5]] = p[0]
                count[p[5]] = [set(), set(), set(), set()]

    seen = set()
    with open(args.i, 'r') as fin:
        for l in fin:
            p=l.strip().split("\t")
            if p[2] not in ['Archaea', 'Bacteria']:
                continue
            seen.add(p[0])
            for i in range(4):
                if len(p) < 6:
                    sys.stderr.write("Not enough elements in {}\n".format("|".join(p)))
                    continue
                count[p[0]][i].add(p[i+2])

    genomeavs = {}
    for i in seen:
        g = genome[i]
        if g not in genomeavs:
            genomeavs[g] = [[], [], [], []]
        for j in range(4):
            genomeavs[g][j].append(len(count[i][j]))
    for g in genomeavs:
        if g not in bodysite and args.b:
            continue
        sys.stdout.write(g)
        if g in bodysite:
            sys.stdout.write("\t{}".format(bodysite[g]))
        else:
            sys.stderr.write("No body site for {}\n".format(g))
            sys.stdout.write("\t")
        for i in range(4):
            av = 1.0 * sum(genomeavs[g][i])/len(genomeavs[g][i])
            sys.stdout.write("\t{}".format(av))
        sys.stdout.write("\n")