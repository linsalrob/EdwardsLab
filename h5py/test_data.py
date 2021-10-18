"""
Create an hdf5 test data set for turbocor


"""

import os
import sys
import argparse
import h5py
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# this is taken from SO: https://stackoverflow.com/questions/18683821/generating-random-correlated-x-and-y-points-using-numpy



xx = np.array([-0.51, 51.2])
yy = np.array([0.33, 51.6])
means = [xx.mean(), yy.mean()]
stds = [xx.std() / 3, yy.std() / 3]
corr = 0.9         # correlation
covs = [[stds[0]**2          , stds[0]*stds[1]*corr],
        [stds[0]*stds[1]*corr,           stds[1]**2]]

m = np.random.multivariate_normal(means, covs, 1000).T
with h5py.File('correlated.h5', 'a') as f:
    if 'data' in f:
        d = f['data']
        d.resize(d.shape[0]+2, axis=0)
        d[-2:] = m
    else:
        f.create_dataset("data", data=m, maxshape=(None,1000), chunks=True)

# sns.scatterplot(m[0], m[1])
# plt.show()

