import matplotlib.pyplot as plt
from scipy.stats.kde import gaussian_kde
import numpy as np
from sklearn.grid_search import GridSearchCV
from sklearn.kernel_ridge import KernelRidge
import sys

x = []
y = []
with open('/home/redwards/Desktop/genus_species_analysis/ecoli_coverage.tsv', 'r') as fin:
    for l in fin:
        p = l.strip().split("\t")
        x.append(float(p[0]))
        y.append(float(p[1]))

px = []
py = []
with open('/home/redwards/Desktop/genus_species_analysis/pseudo_coverage.txt', 'r') as fin:
    for l in fin:
        p = l.strip().split("\t")
        px.append(float(p[0]))
        py.append(float(p[1]))

ny = np.array(y)
nx = np.array(x)
pnx = np.array(px)
pny = np.array(py)


kr = KernelRidge(kernel='rbf', gamma=7.5e-5, alpha=0.001)
kr.fit(nx[:, None], ny[:, None])

x_pred = np.linspace(min(x), max(x), 10000)[:, None]
y_pred = kr.predict(x_pred)


kr.fit(pnx[:, None], pny[:, None])
px_pred = np.linspace(min(px), max(px), 10000)[:, None]
py_pred = kr.predict(px_pred)

fig = plt.figure()
ax = fig.add_subplot(111)


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
    [66, 99], [137, 242],
    [433, 497], [576, 682],
    [822, 879], [986, 1043],
    [1117, 1173], [1243, 1294]
]

illumina = [
    [517, 809],
]


maxy = max([max(y), max(py), max(y_pred), max(py_pred)])

for r in illumina:
    for ix in range(r[0], r[1]):
        ax.bar(ix, maxy, color='lightgrey', edgecolor='lightgrey')

for r in regions:
    for ix in range(r[0], r[1]):
        ax.bar(ix, maxy, color='lightblue', edgecolor='lightblue')

#ax.plot(x, y, color='r')

# ax.plot(px_pred, py_pred, color='red', label='Pseudomonas')
# ax.plot(x_pred, y_pred, color='blue', label='E. coli')
ax.plot(px, py, color='red', label='Pseudomonas')
ax.plot(x, y, color='blue', label='E. coli')
ax.set_ylim([0, maxy])

ax.set_xlabel("Position in the E. coli 16S gene")
ax.set_ylabel("Coverage")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.legend()
fig.set_facecolor('white')

plt.show()
# plt.savefig(args.o)
