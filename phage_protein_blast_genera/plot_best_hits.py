"""
Plot the best hits and their taxonomy

"""

import os
import sys
import argparse
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

xcol = 2 # all proteins: 1 proteins with hit: 2
ycol = 7 # Kingdom : 5 Phylum : 6  Genus : 7   Species : 8
ytitle = 'Genera'


alllocations = ['Appendix', 'Ear', 'Gut', 'Lungs', 'Mouth', 'Nose', 'Skin', 'Vagina']
wantedlocations = ['Gut', 'Mouth']
#wantedlocations = alllocations

skip = ['Escherichia coli']

x = {w:[] for w in wantedlocations}
y = {w:[] for w in wantedlocations}
x['other'] = []
y['other'] = []

with open('/home/redwards/Desktop/genome_best_hits.nums.txt', 'r') as fin:
    for l in fin:
        if l.startswith('#'):
            continue
        p=l.strip().split("\t")
        if p[3] in skip:
            continue
        if p[4] in wantedlocations:
            x[p[4]].append(float(p[xcol]))
            y[p[4]].append(float(p[ycol]))
        else:
            x['other'].append(float(p[xcol]))
            y['other'].append(float(p[ycol]))

fig = plt.figure()
ax = fig.add_subplot(111)

legend = []

ax.plot(x['other'], y['other'], color='0.75', marker='o', linestyle='None')
legend = ['All']

ax.set_xlabel("Number of proteins with similarity")
ax.set_ylabel("Number of different {}".format(ytitle))


colors = ['lightskyblue', 'mediumaquamarine', 'blue', 'green', 'red', 'cyan', 'magenta', 'coral', 'gold', 'aqua']
for i,l in enumerate(wantedlocations):
    dots = ax.plot(x[l], y[l], color=colors[i], marker='o', linestyle='None')
    legend.append(l)

# ax.legend(legend, wantedlocations, numpoints=1)
ax.legend(legend)


# ax.set_xticks(ticks)
# ax.set_xticklabels(labels, rotation='vertical')
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
fig.set_facecolor('white')

plt.tight_layout()
plt.show()
#fig.savefig(figf)


