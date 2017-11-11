"""
Count the normalized bit score for ORF calls

"""

import os
import sys
import argparse
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

def read_ids(idf):
    """
    Read the ID file
    :param idf: an id file that has HOW\tGene
    :return: a hash of gene->how called
    """

    ids={}
    with open(idf, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            if p[1] in ids:
                sys.stderr.write('Error: found {} more than once in {}\n'.format(p[1], idf))
            ids[p[1]]=p[0]
    return ids

def normalized_bitscore(m8file):
    """
    Parse a rapsearch m8 file and calculate the normalized bit score for each entry
    :param m8file:
    :return: a hash of subject id->normalized bit score
    """

    nbs = {}
    with open(m8file, 'r') as f:
        for l in f:
            if l.startswith("#"):
                continue
            # Fields: Query Subject identity        aln-len mismatch        gap-openings    q.start q.end   s.start s.end   log(e-value)    bit-score
            p = l.strip().split("\t")
            if float(p[10]) > -5:
                continue
            s = p[1]
            l = int(p[3])
            b = float(p[11])
            nbs[p[1]] = b/l
            # nbs[p[1]] = b

    return nbs


def read_directory(dir):
    """
    Read a directory of files and put all the data in memory (eek!)
    :param dir: Directory of files to read
    :return:
    """

    data = {}
    for d in os.listdir(dir):
        nbs = normalized_bitscore(os.path.join(dir, d))
        data.update(nbs)
    return data


def separate(nbs, ids):
    """
    Separate the normalized bitscore into all predictions, thea predictions, no predictions
    :param nbs: hash of normalized bit score
    :param ids: hash of IDs
    :return:
    """

    allp = []
    thea = []
    nonp = []

    for i in nbs:
        if i not in ids:
            continue
        if ids[i] == "ANY":
            allp.append(nbs[i])
        elif ids[i] == "THEA":
            thea.append(nbs[i])
        elif ids[i] == "NONE":
            nonp.append(nbs[i])
        else:
            sys.stderr.write("Not really sure what {} is for {}\n".format(ids[i], i))

    return allp, thea, nonp

def plot_nbs(allp, thea, nop, figf):
    """
    Plot the data for all predictions, theapredictions and nopredictions
    :param allp: all predictions
    :param thea: thea predictions
    :param nop: orf predictions (not called)
    :param figf: the output file name for the figure
    :return:
    """

    alldata = [allp, thea, nop]
    labels = ["All", "THEA", "None"]

    sys.stderr.write("Lengths: All: {} THEA: {} NONE: {}\n".format(len(allp), len(thea), len(nop)))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    #ax.boxplot(alldata)
    ax.violinplot(alldata, [1,2,3], showmeans=True)
    ax.set_xlabel("Predictions")
    ax.set_ylabel("Normalized bit score")
    ax.set_xticks([1,2,3])
    ax.set_xticklabels(labels)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    fig.set_facecolor('white')

    # plt.show()
    fig.savefig(figf)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate the normalized bit score")
    parser.add_argument('-m', help='m8 output file from rapsearch')
    parser.add_argument('-d', help='directory of m8 files')
    parser.add_argument('-i', help='ids file', required=True)
    parser.add_argument('-f', help='figure output file name', default="fig.png")
    parser.add_argument('-v', help='verbose output')
    args = parser.parse_args()

    ids = read_ids(args.i)
    nbs = {}
    if args.m:
        nbs = normalized_bitscore(args.m)
    elif args.d:
        nbs = read_directory(args.d)
    else:
        sys.stderr.write("ERROR: either -d or -m must be supplied. Use -h for help")
        sys.exit(-1)
    a, t, n = separate(nbs, ids)
    plot_nbs(a,t,n,args.f)
