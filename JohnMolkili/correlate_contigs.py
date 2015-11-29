import os
import sys

from scipy.stats.stats import pearsonr
import argparse

__author__ = 'Rob Edwards'


def merge_clust(c1, c2, inclust, clustermembers):
    if c2 < c1:
        [c1, c2] = [c2, c1]

    for c in clustermembers[c2]:
        clustermembers[c1].add(c)
        inclust[c] = c1

    clustermembers.pop(c2)
    return inclust, clustermembers

def notzero(x): return x > 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='calculate pairwise pearson correlations between contigs and then cluster them')
    parser.add_argument('-d', help='data table with contigs in rows and occurence in columns', required=True)
    parser.add_argument('-t', help='pearson correlation minimum for clustering (default = 0.9)', default=0.9, type=float)
    args = parser.parse_args()

    data = {}
    headers = []
    with open(args.d, 'r') as f:
        for l in f:
            p=l.strip().split("\t")
            if headers == []:
                headers = p
            else:
                data[p[0]] = map(int, p[1:])

    allcontigs = data.keys()
    allcontigs.sort()
    dist = {x: {} for x in allcontigs}

    # we're just going to use contigs where at least 10 samples are not zero
    nonzero = []
    for c in allcontigs:
        nz = filter(notzero, data[c])
        if len(nz) > 10:
            nonzero.append(c)

    sys.stderr.write("Before filtering we had {} contigs, after filtering on 10 > 0 we have {} contigs\n".format(len(allcontigs), len(nonzero)))

    cluster = 0
    inclust = {}
    clustermembers = {}

    for i in range(len(nonzero)):
        cfr = nonzero[i]
        for j in range(i, len(nonzero)):
            cto = nonzero[j]
            if cfr in inclust and cto in inclust and inclust[cfr] == inclust[cto]:
                # no need to calculate!
                continue
            dist[cfr][cto] = dist[cto][cfr] = pearsonr(data[cfr], data[cto])
            if dist[cfr][cto] < args.t:
                if cfr in inclust and cto in inclust:
                    # need to merge these two clusters
                    inclust, clustermembers = merge_clust(inclust[cfr], inclust[cto], inclust, clustermembers)
                elif cfr in inclust:
                    inclust[cto] = inclust[cfr]
                    clustermembers[inclust[cfr]].append(cto)
                elif cto in inclust:
                    inclust[cfr] = inclust[cto]
                    clustermembers[inclust[cto]].append(cfr)
                else:
                    inclust[cfr] = cluster
                    inclust[cto] = cluster
                    clustermembers[cluster] = {cfr, cto}
                    cluster += 1


    for contig in inclust:
        print("{}\t{}".format(contig, inclust[contig]))
