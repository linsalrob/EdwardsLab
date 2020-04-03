"""
Count adjacent orfs in the blast output files.
"""

import os
import sys
import argparse


def count_adjacent_orfs(sample, fastafile, blastfile, adjacentout, nohitsout, searchtype):
    """
    Count hits where two adjacent orfs match
    to the same protein
    """
    sys.stderr.write(f"{bcolors.GREEN}Counting adjacent ORFs for {sample} and {searchtype}{bcolors.ENDC}\n")
    orfs = []
    with open(fastafile, 'r') as f:
        for l in f:
            if l.startswith('>'):
                p = l.split()
                seqid = p[0].replace('>', '', 1)
                orfs.append(seqid)
    hits = {}
    with open(blastfile, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            if p[0] not in hits:
                hits[p[0]] = set()
            hits[p[0]].add(p[1])

    adjacent = 0
    nohits = 0
    for i, o in enumerate(orfs):
        if i >= len(orfs) - 1:
            continue
        if o not in hits:
            nohits += 1
            continue
        if orfs[i + 1] not in hits:
            continue
        for h in hits[o]:
            if h in hits[orfs[i + 1]]:
                adjacent += 1
                break
    with open(adjacentout, 'w') as out:
        out.write(f"{sample}\tNumber of orfs with adjacent {searchtype} similarities\t")
        out.write(f"[adjacent orfs, total orfs, fraction adjacent]\t")
        out.write(f"{adjacent}\t{len(orfs)}\t{adjacent / len(orfs)}\n")
    with open(nohitsout, 'w') as out:
        out.write(f"{sample}\tNumber of orfs with no hits\t")
        out.write(f"[nohits, total orfs, fraction no hits]\t")
        out.write(f"{nohits}\t{len(orfs)}\t{nohits / len(orfs)}\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-s', help='sample name used in output', required=True)
    parser.add_argument('-f', help='fasta file of protein sequences', required=True)
    parser.add_argument('-b', help='blast m8 file', required=True)
    parser.add_argument('-a', help='adjacent orfs output file', required=True)
    parser.add_argument('-n', help='no hits output file', required=True)
    parser.add_argument('-t', help='search type  (e.g. phage, bacteria) (used in output)', required=True)
    args = parser.parse_args()

    count_adjacent_orfs(args.s, args.f, args.b, args.a, args.n, args.t)