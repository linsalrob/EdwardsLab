"""
Plot a 3D scatter plot of the alignment scores from alignment_score.py
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

import tkinter as tk
from tkinter import filedialog

if __name__ == '__main__':
    """
    parser = argparse.ArgumentParser(description="Plot a 3D scatter plot of the alignment scores from alignment_score.py")
    parser.add_argument('-f', help='text separated data file output from alignment_score.py', required=True)
    args = parser.parse_args()

    filename = args.f
    """

    root = tk.Tk()
    root.withdraw()
    filename = filedialog.askopenfilename()

    x = []
    y = []
    z = []
    legends = ['seq1', 'seq2', '1-mer', '2-mer', '3-mer']
    with open(filename, 'r') as fin:
        for l in fin:
            p = l.strip().split("\t")
            x.append(float(p[2]))
            y.append(float(p[3]))
            z.append(float(p[4]))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(x, y, z)
    ax.legend()

    ax.set_xlim(0, max(x))
    ax.set_ylim(0, max(y))
    ax.set_zlim(0, max(z))
    ax.set_xlabel(legends[2])
    ax.set_ylabel(legends[3])
    ax.set_zlabel(legends[4])
    pickle.dump(fig, open('/home/redwards/Desktop/3dfig.pickle', 'wb'))
    plt.show()

