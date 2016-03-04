"""
Plot the SRA data from PARTIE. Espeically plot the percent phage + prok against percent 16S
"""

import os
import random
import sys

import matplotlib.lines
import argparse
import matplotlib.colors
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot the SRA data from PARTIE")
    parser.add_argument('-p', help='partie output file', required=True)
    parser.add_argument('-e', help='Experiment library file (run id\texperiment library', required=True)
    args = parser.parse_args()

    explib = {}
    with open(args.e, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            if len(p) == 1:
                p.append('OTHER')
            explib[p[0]] = p[1]

    data = {}
    experimentlibraries = {}
    with open(args.p, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            # data is unique kmers, percent 16S, percent phage, percent prok, percent prok + phage
            p[0] = p[0].replace('.sra', '')
            data[p[0]] = [float(p[2]), float(p[4]), float(p[6]), float(p[8]), float(p[6]) + float(p[8])]
            if p[0] in explib:
                if explib[p[0]] in experimentlibraries:
                    experimentlibraries[explib[p[0]]].append(p[0])
                else:
                    experimentlibraries[explib[p[0]]] = [p[0]]
            else:
                sys.stderr.write("No " + p[0] + " in exp\n")

    sra_ids = data.keys()

    # allcolors = matplotlib.colors.cnames.keys()
    allcolors = ['indigo', 'gold', 'hotpink', 'firebrick', 'indianred', 'yellow',
                 'mistyrose', 'olive', 'pink', 'tomato', 'orangered', 'navajowhite', 'lime', 'palegreen', 'greenyellow',
                 'burlywood', 'seashell', 'mediumspringgreen', 'fuchsia', 'papayawhip', 'blanchedalmond', 'chartreuse',
                 'dimgray', 'black', 'peachpuff', 'springgreen', 'aquamarine', 'white', 'orange', 'brown', 'ivory',
                 'dodgerblue', 'peru', 'lawngreen', 'chocolate', 'crimson', 'forestgreen', 'slateblue', 'cyan',
                 'mintcream', 'silver', 'antiquewhite', 'mediumorchid', 'skyblue', 'gray', 'goldenrod', 'floralwhite',
                 'moccasin', 'saddlebrown', 'grey', 'mediumvioletred', 'slategrey', 'red', 'deeppink', 'limegreen',
                 'palegoldenrod', 'plum', 'turquoise', 'lavender', 'maroon', 'yellowgreen', 'sandybrown', 'thistle',
                 'violet', 'navy', 'magenta', 'dimgrey', 'tan', 'rosybrown', 'olivedrab', 'blue', 'ghostwhite',
                 'honeydew', 'cornflowerblue', 'linen', 'powderblue', 'seagreen', 'snow', 'sienna', 'mediumblue',
                 'royalblue', 'green', 'mediumpurple', 'midnightblue', 'cornsilk', 'paleturquoise', 'bisque',
                 'slategray', 'khaki', 'wheat', 'teal', 'deepskyblue', 'salmon', 'steelblue', 'palevioletred',
                 'aliceblue', 'orchid', 'gainsboro', 'mediumseagreen', 'mediumturquoise', 'lemonchiffon', 'cadetblue',
                 'lavenderblush', 'coral', 'purple', 'aqua', 'whitesmoke', 'mediumslateblue', 'mediumaquamarine',
                 'beige', 'blueviolet', 'azure', 'oldlace']


    labels = {}

    phageprok = []
    sixs = []
    kms = []

    #for e in experimentlibraries:
    libs2use = ['WGS', 'AMPLICON', 'CLONE', 'OTHER', 'RNA-Seq', 'WGA']
    for e in libs2use:
        col = allcolors.pop(0)
        print(e + "\t" + str(len(experimentlibraries[e])))


        labels[col] = e
        phageprokn = []
        sixsn = []
        kmsn = []
        """
        for i in range(1000):
            r = random.randint(0, len(experimentlibraries[e]) - 1)
            phageprokn.append(data[experimentlibraries[e][r]][3])
            sixsn.append(data[experimentlibraries[e][r]][1])
            kmsn.append(data[experimentlibraries[e][r]][0])
        """
        for r in range(len(experimentlibraries[e])):
            phageprokn.append(data[experimentlibraries[e][r]][3])
            sixsn.append(data[experimentlibraries[e][r]][1])
            kmsn.append(data[experimentlibraries[e][r]][0])

        phageprok.append(phageprokn)
        sixs.append(sixsn)
        kms.append(kmsn)


    f, (ppax, sxax, kmax) = plt.subplots(1,3)

    ppax.boxplot(phageprok)
    plt.setp(ppax, xticks=[1, 2, 3, 4, 5, 6], xticklabels=libs2use)

    sxax.boxplot(sixs)
    plt.setp(sxax, xticks=[1, 2, 3, 4, 5, 6], xticklabels=libs2use)

    kmax.boxplot(kms)
    plt.setp(kmax, xticks=[1, 2, 3, 4, 5, 6], xticklabels=libs2use)

    ppax.set_ylabel("Percent Prokaryote")
    sxax.set_ylabel("Percent 16S")
    kmax.set_ylabel("KMERs")

    # this is to get the legend on a 3D plot
    scatterproxy = []
    labeltexts = []
    for color in labels:
        scatterproxy.append(matplotlib.lines.Line2D([0], [0], linestyle="none", c=color, marker='o'))
        labeltexts.append(labels[color])


    plt.show()
