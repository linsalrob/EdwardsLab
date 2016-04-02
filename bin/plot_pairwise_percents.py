"""
Plot the pairwise percents.

The perl code perl average_pairwise.pl generates a JSON data structure with the percent IDs from the same organisms at the different levels

This plots that

"""
import json

import matplotlib.pyplot as plt


def median(lst):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = (lstLen - 1) // 2 # // is the floor division

    if (lstLen % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1]) / 2.0



filename = '/home/redwards/Desktop/identical_percent_ids.tsv'
with open(filename, 'r') as f:
    data = json.load(f)


tax = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
alldata = []
for t in tax:
    # data[t] = map(float, data[t])
    floatdata = map(float, data[t])
    alldata.append(floatdata)
    print("{}\t{}\t{}\t{}".format(t, len(floatdata), 1.0*sum(floatdata)/len(floatdata), median(floatdata)))


fig = plt.figure()
ax = fig.add_subplot(111)

ax.boxplot(alldata)
ax.set_xlabel("Phylogeny")
ax.set_ylabel("Average percent identity")
ax.set_xticklabels(tax)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
fig.set_facecolor('white')

plt.show()

