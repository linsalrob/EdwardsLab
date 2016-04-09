import matplotlib.pyplot as plt
from scipy.stats.kde import gaussian_kde
import numpy as np
from sklearn.grid_search import GridSearchCV
from sklearn.neighbors import KernelDensity
import sys

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

ny = np.array(y)
nx = np.array(x)

xbins = []
ybins = []

for i in range(0, len(y), 20):
    xbins.append(i)
    ybins.append(0)
    for j in range(i, i+19):
        if j >= len(y):
            break
        ybins[-1]+=y[j]
    ybins[-1] = ybins[-1]/20

nybins = np.array(ybins)

bandwidth = 0.2
if 0:
    kde_skl = KernelDensity(bandwidth=bandwidth)
    kde_skl.fit(nybins[:, np.newaxis])

    # score_samples() returns the log-likelihood of the samples
    log_pdf = kde_skl.score_samples(nx[:, np.newaxis])
    density = np.exp(log_pdf)


grid = GridSearchCV(KernelDensity(kernel='gaussian'),
                    {'bandwidth': np.linspace(0.1, 2000, 200)},
                    cv=31)  # 20-fold cross-validation

grid.fit(nybins[:, np.newaxis])
print(grid.best_params_)

kde = grid.best_estimator_
pdf = np.exp(kde.score_samples(ny[:, None]))


print("X: " + str(x[0:10]))


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
    for r in illumina:
        for ix in range(r[0], r[1]):
            ax.bar(ix, maxy, color='lightgrey', edgecolor='lightgrey')

    for r in regions:
        for ix in range(r[0], r[1]):
            ax.bar(ix, maxy, color='lightblue', edgecolor='lightblue')

ax.plot(x, y, color='r')

ax.bar(xbins, ybins, width=20, alpha=0.5, color='lightgreen')

ax2.plot(x, pdf, color='blue')
ax2.set_ylabel("likelihood score")

ax.set_xlabel("Position in the E. coli 16S gene")
ax.set_ylabel("Coverage")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

fig.set_facecolor('white')

plt.show()
# plt.savefig(args.o)
