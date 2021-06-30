"""
We have created mash sketches of the GPDB database, the MGV database, and the SDSU phage, and
this will figure out the top hits and summarize their familes.
"""

import os
import sys
import argparse


def best_hits(distf, maxscore, verbose=False):
    """
    Find the best hits
    """
    bh = {}
    allph = set()
    with open(distf, 'r') as din:
        for li in din:
            p = li.strip().split("\t")
            if float(p[3]) <= maxscore:
                if p[0] not in bh:
                    bh[p[0]] = set()
                bh[p[0]].add(p[1])
            allph.add(p[0])

    if verbose:
        for p in allph:
            if p not in bh:
                sys.stderr.write(f"WARNING: With a score of {maxscore} did not find any hits to {p}\n")
    return bh

def find_vc(mdf, genomecol, vccol, verbose=False):
    """
    Read the metadata file and return a hash of genome->viral cluster
    """
    vc = {}
    with open(mdf, 'r') as fin:
        for li in fin:
            p = li.strip().split("\t")
            vc[p[genomecol]] = p[vccol]
    if verbose:
        sys.stderr.write(f"Found {len(vc)} virus clusters in {mdf}\n")
    return vc


def count_hits(bh, vc, verbose=False):
    """
    Count the vc hits per genome
    """

    hc = {}
    for g in bh:
        hc[g] = {}
        for b in bh[g]:
            hc[g][vc[b]] = hc[g].get(vc[b], 0) + 1
        besthit = None
        bhc = 0
        for h in hc[g]:
            if hc[g][h] > bhc:
                bhc = hc[g][h]
                besthit = h
        print(f"{g}\t{besthit}\t{bhc}\t{len(bh[g])}")

    return hc


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-d', help='mash distance file', required=True)
    parser.add_argument('-c', help='distance cutoff score, default = 0', default=0, type=float)
    parser.add_argument('-m', help='metadata file', required=True)
    parser.add_argument('-g', help='genome column, default = 0', default=0, type=int)
    parser.add_argument('-l', help='virus cluster col in the metadata file', type=int, required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    bh = best_hits(args.d, args.c, args.v)
    vc = find_vc(args.m, args.g, args.l, args.v)
    count_hits(bh, vc,args.v)
