import matplotlib.pyplot as plt
from scipy.stats.kde import gaussian_kde
import numpy as np
from sklearn.grid_search import GridSearchCV
from sklearn.kernel_ridge import KernelRidge
import sys

x = []
y = []
with open('/home/redwards/Desktop/genus_species_analysis/ecoli_coverage.tsv', 'r') as fin:
#with open('/home/redwards/Desktop/genus_species_analysis/pseudo_coverage.txt', 'r') as fin:
    for l in fin:
      p=l.strip().split("\t")
      x.append(float(p[0]))
      y.append(float(p[1]))


ny = np.array(y)
nx = np.array(x)

grid = GridSearchCV(
    KernelRidge(kernel='rbf', gamma=1e-4),
    param_grid={"alpha": [0.1, 0.01, 0.001]},
                cv=5)  # 20-fold cross-validation

# param_grid={"alpha": np.logspace(-10, 10, 10),
# "gamma": np.logspace(-4, -3, 5)},

grid.fit(nx[:, None], ny[:, None])
print(grid.best_params_)

xaxis_pred = np.linspace(min(x), max(x), 10000)[:, None]
yaxis_pred = grid.predict(xaxis_pred)



fig = plt.figure()
ax = fig.add_subplot(111)

ax2=ax.twinx()


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

if 0:
    maxy = max(y)
    for r in illumina:
        for ix in range(r[0], r[1]):
            ax.bar(ix, maxy, color='lightgrey', edgecolor='lightgrey')

    for r in regions:
        for ix in range(r[0], r[1]):
            ax.bar(ix, maxy, color='lightblue', edgecolor='lightblue')

ax.plot(x, y, color='r')

ax2.plot(xaxis_pred, yaxis_pred, color='blue')
ax2.set_ylabel("predictions")
ax2.set_ylim([0, 5000])

ax.set_xlabel("Position in the E. coli 16S gene")
ax.set_ylabel("Coverage")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

fig.set_facecolor('white')

plt.show()
# plt.savefig(args.o)
