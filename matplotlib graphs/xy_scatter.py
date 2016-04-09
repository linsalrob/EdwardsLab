"""
Create an x-y scatter plot and save it as an SVG
"""

import os
import sys

import argparse

import matplotlib

#matplotlib.use('SVG')
import matplotlib.pyplot as plt

def read_file(filename, first_line_labels):
    """
    Read the tsv file and return the labels and the data.
    :param filename: The data file
    :type filename: str
    :param labels: Whether the first line has labels
    :type labels: bool
    :return: array and 2-D array of labels and data
    :rtype: array, array
    """

    labels = ['X', 'Y']

    x=[]
    y=[]
    with open(filename, 'r') as fin:
        if first_line_labels:
            labels = fin.readline().strip().split("\t")

        for l in fin:
            p=l.strip().split("\t")
            x.append(float(p[0]))
            y.append(float(p[1]))

    return labels, [x, y]


if __name__ == '__main__':
    if (1):
        parser = argparse.ArgumentParser(description="Generate an x-y scatter plot from a tsv file")
        parser.add_argument('-f', help='tsv file with data', required=True)
        parser.add_argument('-o', help="output file name", required=True)
        parser.add_argument('-l', help='first line of tsv file has labels', action='store_true')
        parser.add_argument('-c', help='line color. See http://matplotlib.org/api/colors_api.html for codes. Default=black', default='k')
        parser.add_argument('-b', help='add break lines that go from 0 to max(y) at position(s) provided', action='append')
        args = parser.parse_args()

        labels, data = read_file(args.f, args.l)

    # labels, data = read_file('/home/redwards/Desktop/integrase_clusters.tsv', True)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    maxy=max(data[1])
    #for l in [10,20,50,75]:
    for l in args.b:
        x = [l, l]
        y = [0, maxy]
        ax.plot(x, y, color='lightgrey')

    ax.plot(data[0], data[1], color=args.c)
    # ax.plot(data[0], data[1], color='k')

    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()


    fig.set_facecolor('white')

    #plt.show()
    plt.savefig(args.o)

