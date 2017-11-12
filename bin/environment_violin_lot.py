"""

"""

import os
import sys
import argparse
import gzip

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot(f, figf):
    freq = {}
    with open(f, 'r') as i:
        for l in i:
            p=l.strip().split("\t")

            if 'human' in p[0]:
                p[0]='human'
            if 'human' in p[1]:
                p[1]='human'
            if p[0] < p[1]:
                (g1, g2) = (p[0], p[1])
            else:
                (g2, g1) = (p[0], p[1])

            if g1 not in freq:
                freq[g1] = {}
            if g2 not in freq[g1]:
                freq[g1][g2] = []

            freq[g1][g2].append(float(p[2]))

    labels = []
    scores = []
    sames = []
    count = 1
    ticks = []
    for g1 in freq.keys():
        for g2 in freq[g1].keys():
            if len(freq[g1][g2]) < 1000000:
                continue
            labels.append("{}-{}".format(g1, g2))
            scores.append(freq[g1][g2])
            if g1 == g2:
                sames.append(1)
            else:
                sames.append(0)
            ticks.append(count)
            count += 1



    fig = plt.figure()
    ax = fig.add_subplot(111)

    # ax.boxplot(alldata)
    vp = ax.violinplot(scores, ticks, showmeans=True)
    for i, j in enumerate(vp['bodies']):
        if sames[i]:
            j.set_color('red')

    ax.set_xlabel("Environments")
    ax.set_ylabel("Genome Distances")
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, rotation='vertical')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    fig.set_facecolor('white')

    plt.tight_layout()
    #plt.show()
    fig.savefig(figf)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Draw a voilin plot of same and different environments")
    parser.add_argument('-f', help='file of environments and distances', default="/home/redwards/Desktop/first_part.tsv")
    parser.add_argument('-o', help='output image file', default='/home/redwards/Desktop/out.png')
    parser.add_argument('-v', help='verbose output')
    args = parser.parse_args()

    plot(args.f, args.o)