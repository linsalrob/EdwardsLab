

import os
import sys
import re
import matplotlib.pyplot as plt
import numpy as np
import roblib


longest_sequence = 1600
min_len = 45
hits = {}

subtract_medians = False
subtract_means = False

min_coverage = 0 # miniumum number of bases that must be in 16S gene to be included

bacteria = {"158878.1": "Staphylococcus aureus subsp. aureus Mu50", "160490.1": "Streptococcus pyogenes M1 GAS",
            "195099.3": "Campylobacter jejuni RM1221", "196164.1": "Corynebacterium efficiens YS-314",
            "208964.1": "Pseudomonas aeruginosa PAO1", "226185.1": "Enterococcus faecalis V583",
            "264731.4": "Prevotella ruminicola 23", "411466.7": "Actinomyces odontolyticus ATCC 17982",
            "479436.6": "Veillonella parvula DSM 2008", "592010.4": "Abiotrophia defectiva ATCC 49176",
            "762948.4": "Rothia dentocariosa ATCC 17931", "83333.1": "Escherichia coli K12"}


def initiate_seq(sid):
    """
    Just add this sequence to the hits if it is not already there

    :param sid: sequence id
    :type sid: str
    """

    global hits

    if sid in hits:
        return

    hits[sid] = []
    for i in range(longest_sequence):
        hits[sid].append(0)


def plot_vregions(ax, maxy):
    """
    Plot the v regions
    :param ax: matplotlib figure
    :type ax: matplotlib figure
    """

    """
    These regions come from http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2562909/
    v1: 66-99
    v2: 137-242
    v3: 433-497
    v4: 576-682
    v5: 822-879
    v6: 986-1043
    v7: 1117-1173
    v8: 1243-1294

    """

    illumina = [
        [517, 809],
    ]

    regions = [
        [66, 99], [137, 242],
        [433, 497], [576, 682],
        [822, 879], [986, 1043],
        [1117, 1173], [1243, 1294]
    ]

    for r in illumina:
        for ix in range(r[0], r[1]):
            ax.bar(ix, maxy, color='lightgrey', edgecolor='lightgrey')

    for r in regions:
        for ix in range(r[0], r[1]):
            ax.bar(ix, maxy, color='lightblue', edgecolor='lightblue')


def notzero(x): return x > 0


def clean_hits():
    """
    Remove hits with < coverage than minimum
    :return:
    :rtype:
    """

    temp = {}
    for s in hits:
        data = filter(notzero, hits[s])
        if len(data) > min_coverage:
            temp[s] = hits[s]
    return temp


#with open('/home/redwards/Desktop/16S_all.cf94.blastn', 'r') as bin:
with open('/home/redwards/Desktop/16S_all_padded.blastn', 'r') as bin:
    for l in bin:
        p=l.strip().split('\t')
        for i in [2, 10, 11]:
            p[i] = float(p[i])
        for i in range(3,10):
            p[i] = int(p[i])

        if p[3] < min_len:
            continue

        if p[8] > p[9]:
            (p[8], p[9]) = (p[9], p[8])

        m=re.match('fig\|(\d+\.\d+)', p[1])
        if not m:
            sys.stderr.write("Can't parse" + p[1] + "\n")
            continue

        sid = bacteria[m.group(1)]
        if sid not in hits:
            initiate_seq(sid)

        for i in range(p[8], p[9]+1):
            hits[sid][i] = hits[sid][i]+1


# flter for coverage
sys.stderr.write("Before cleaning there are " + str(len(hits.keys())) + " hits\n")
hits = clean_hits()
sys.stderr.write("After cleaning there are " + str(len(hits.keys())) + " hits\n")

# calculate the median for each position

sys.stderr.write("calculating medians\n")
median = []
mean = []
for i in range(longest_sequence):
    data = [hits[s][i] for s in hits]
    median.append(roblib.median(data))
    mean.append(roblib.mean(data))

median = np.array(median)
mean = np.array(mean)

if subtract_medians:
    for s in hits:
        hits[s] = np.array(hits[s]) - median

if subtract_means:
    for s in hits:
        hits[s] = np.array(hits[s]) - mean

maxy = 0
for s in hits:
    m = max(hits[s])
    if m > maxy: maxy = m

# now we just need to make a plot for each of the sequences
fig = plt.figure()
ax = fig.add_subplot(111)


plot_vregions(ax, maxy)


x=range(longest_sequence)

for s in hits:
    ax.plot(x, hits[s], label=s)

ax.set_xlabel("Position in the 16S gene")
ax.set_ylabel("Coverage")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

ax.legend()

fig.set_facecolor('white')

plt.show()





