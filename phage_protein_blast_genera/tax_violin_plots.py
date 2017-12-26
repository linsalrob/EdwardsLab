"""

"""

import os
import sys
import argparse

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-f', help='Genome average output file (from genera_per_phage_protein.py', default='/home/redwards/Desktop/gav.out')
    parser.add_argument('-n', help='taxonomy name one of: kingdom / phylum / genus / species', default='phylum')
    parser.add_argument('-v', help='verbose output', action="store_true")

    args = parser.parse_args()

    col = None
    colkey = {'kingdom' : 2, 'phylum' : 3, 'genus' : 4, 'species' : 5}
    if args.n not in colkey:
        sys.stderr.write("Sorry, taxonomy name must be one of {}\n".format("|".join(list(colkey.keys()))))
        sys.exit(-1)
    col = colkey[args.n]

    want = {'Gut', 'Mouth'}

    data = {}
    with open(args.f, 'r') as fin:
        for l in fin:
            p=l.strip().split("\t")
            if p[1] not in want:
                True
                continue
            if p[1] not in data:
                data[p[1]] = []
            data[p[1]].append(float(p[col]))

    labels = sorted(data.keys())
    scores = []
    count = 1
    ticks = []
    for l in labels:
        scores.append(data[l])
        ticks.append(count)
        count += 1

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # ax.boxplot(alldata)
    vp = ax.violinplot(scores, showmeans=True)
    ax.set_xlabel("Body Site")
    ax.set_ylabel("Average number of {}".format(args.n))
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, rotation='vertical')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    fig.set_facecolor('white')

    plt.tight_layout()
    plt.show()
    #fig.savefig(figf)
