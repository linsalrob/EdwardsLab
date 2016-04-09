"""
Plot a 3D scatter plot
"""

import os
import random
import sys

import matplotlib.lines
import argparse
import matplotlib.colors
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot a 3D scatter plot")
    parser.add_argument('-f', help='text separated data file', required=True)
    parser.add_argument('-l', help='File has a header line with axis titles', action='store_true')
    args = parser.parse_args()

    titles = []
    x = []
    y = []
    z = []
    legends = ['what', 'x', 'y', 'z']
    with open(args.f, 'r') as fin:
        if args.l:
            legends = fin.readline().strip().split("\t")
        for l in fin:
            p = l.strip().split("\t")
            titles.append(p[0])
            x.append(float(p[1]))
            y.append(float(p[2]))
            z.append(float(p[3]))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(x, y, z)
    ax.legend()

    ax.set_xlim(0, max(x))
    ax.set_ylim(0, max(y))
    ax.set_zlim(0, max(z))
    ax.set_xlabel(legends[1])
    ax.set_ylabel(legends[2])
    ax.set_zlabel(legends[3])
    pickle.dump(fig, open('/home/redwards/Desktop/3dfig.pickle', 'wb'))
    plt.show()

