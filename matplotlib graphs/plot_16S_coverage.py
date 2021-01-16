import matplotlib.pyplot as plt
import numpy

x = []
y = []
with open('/home/redwards/Desktop/genus_species_analysis/ecoli_coverage.tsv', 'r') as fin:
#with open('/home/redwards/Desktop/genus_species_analysis/pseudo_coverage.txt', 'r') as fin:
    for l in fin:
      p=l.strip().split("\t")
      x.append(float(p[0]))
      y.append(float(p[1]))

fig = plt.figure()
ax = fig.add_subplot(111)
maxy = max(y)

ax.plot(x, y, color='r')

ax.plot(xs, kdepdf(xs), color='blue')

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


regions = [
    [66,99], [137, 242],
    [433, 497], [576, 682],
    [822, 879], [986, 1043],
    [1117, 1173], [1243, 1294]
]

illumina = [
    [517, 809],
]

for r in illumina:
    for x in range(r[0], r[1]):
        ax.bar(x, maxy, color='lightgrey', edgecolor='lightgrey')

for r in regions:
    for x in range(r[0], r[1]):
        ax.bar(x, maxy, color='lightblue', edgecolor='lightblue')


ax.set_xlabel("Position in the E. coli 16S gene")
ax.set_ylabel("Coverage")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

fig.set_facecolor('white')

plt.show()