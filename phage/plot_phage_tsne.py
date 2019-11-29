"""
Read the phage fasta and plot a tSNE. This is adapted from phage_protein_data.py
"""

import os
import sys
import argparse
import json

import numpy as np
from sklearn.mixture import GaussianMixture
from sklearn.manifold import TSNE

from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d

import seaborn as sns
sns.set()

import pandas as pd

import itertools

from roblib import bcolors, stream_fasta

__author__ = 'Rob Edwards'


def count_pairwise(faf, kmer, verbose=True):
    """
    Count all pairwise amino acids
    :param faf: fasta file
    :param kmer: kmer size
    :param verbose: more output
    """

    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Reading {faf}{bcolors.ENDC}\n")

    count = {}
    fns = {}
    for sidf, seq in stream_fasta(faf, whole_id=True):
        posn=0
        # split the sequence id into sequence and function
        # may need to provide an alternate way to do this
        # note we also strip out [organism name]
        sid = sidf[:sidf.index(" ")]
        try:
            fns[sid] = sidf[sidf.index(" ")+1:sidf.index("[")-1]
        except:
            fns[sid] = sidf[sidf.index(" ")+1:]

        count[sid]={}
        while posn < len(seq)-(kmer-1):
            count[sid][seq[posn:posn+kmer]] = count[sid].get(seq[posn:posn+kmer], 0)+1
            posn += 1
        # normalize by protein length
        for aa in count[sid]:
            count[sid][aa] /= len(seq)
    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Read {faf}{bcolors.ENDC}\n")

    return count, fns


def count_pairwise_no_fn(faf, kmer, verbose=True):
    """
    Count all pairwise amino acids
    :param faf: fasta file
    :param kmer: kmer size
    :param verbose: more output
    """

    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Reading {faf}{bcolors.ENDC}\n")

    count = {}
    for sid, seq in stream_fasta(faf, whole_id=True):
        posn=0
        count[sid]={}
        while posn < len(seq)-(kmer-1):
            count[sid][seq[posn:posn+kmer]] = count[sid].get(seq[posn:posn+kmer], 0)+1
            posn += 1
        # normalize by protein length
        for aa in count[sid]:
            count[sid][aa] /= len(seq)
    return count


def read_fns(fnf, verbose=True):
    """
    Read the tsv of protein id/function. We use this when the function
    is not included in the fasta file
    """
    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Reading functions from {fnf}{bcolors.ENDC}\n")
    
    fns={}
    with open(fnf, 'r') as f:
        for l in f:
            p=l.strip().split("\t")
            fns[p[0]]=p[1]
    return fns


def counts_to_df(counts, kmer, verbose=False):
    """
    Convert the count dict to a pandas data frame
    :param counts: the dict of sequences, aa and counts
    :param kmer: the kmer size
    :param verbose: more output
    :return pd.DataFrame
    """

    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Making possibilities{bcolors.ENDC}\n")

    # this works for K< 7
    # iupac= ["".join(tple) for tple in itertools.product('ACDEFGHIKLMNPQRSTVWY', repeat=kmer)]

    allmers = set()
    for s in counts:
        allmers.update(counts[s].keys())

    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Converting to matrix {bcolors.ENDC}\n")

    aml = list(allmers)
    aml.sort()

    pc = {}
    for sid in counts:
        pc[sid] = []
        for aa in aml:
            pc[sid].append(counts[sid].get(aa, 0))


    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Converting to dataframe{bcolors.ENDC}\n")

    df = pd.DataFrame.from_dict(pc, orient='index', columns=aml)
    return df

def unique_functions(fns, verbose=False):
    """
    Returns the number of unique functions
    :param functions: the dict of functions
    :param verbose: more output
    :returns: int
    """

    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Unique functions{bcolors.ENDC}\n")

    return len(set(fns.values()))

def fns_to_list(fns, df, verbose=False):
    """
    Conver the functions to an ordered list based on the df
    """

    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Listifying functions{bcolors.ENDC}\n")

    return [fns[x] for x in df.index.values]

def create_tSNE(df, fnar, tout, verbose):
    """
    Generate the tSNE data and write it to a text file
    :param df: the data frame
    :param fnar: the function array
    :param tout: the output file to write
    :param verbose: more output
    :return: the tSNE array
    """

    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Generating tSNE {bcolors.ENDC}\n")

    tsne = TSNE(n_components=3)
    X = tsne.fit_transform(df)
    if tout:
        data = {
            'tSNE' : X.tolist(),
            'functions' : fnar
        }
        with open(tout, 'w') as out:
            json.dump(data, out)
    
    return X


def filternz(df, verbose=True):
    """
    Filter columns that sum to 0. We want to remove these
    :param df: the data frame
    :param verbose: more output
    """
    
    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Filtering to remove column sums == 0{bcolors.ENDC}\n")
    return df.loc[:, df.sum() != 0]

def plot_2d(tsne, fnar, palette, outpng, verbose):
    """
    Create the 2d plot
    :param tsne: tSNE array
    :param fnar: functions list
    :param outpng: base name for pg output
    :return: nothing
    """
    
    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Plotting 2D tSNE{bcolors.ENDC}\n")

    snsplot = sns.scatterplot(tsne[:,0], tsne[:,1], legend="full", hue=fnar, palette=palette)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(outpng + ".png")

def plot_2d_sz(tsne, fnar, palette, outpng, verbose):
    """
    Create the 2d plot
    :param tsne: tSNE array
    :param fnar: functions list
    :param outpng: base name for pg output
    :return: nothing
    """
    
    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Plotting 2D tSNE by size{bcolors.ENDC}\n")

    sp = sns.scatterplot(x=tsne[:,0], y=tsne[:,1], s=tsne[:,2], legend="full", hue=fnar, palette=palette)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(outpng + ".sz.png")

def plot_3d(tsne, fnar, palette, outpng, verbose):
    """
    Create the 3d plot
    :param tsne: tSNE array
    :param fnar: functions list
    :param outpng: base name for pg output
    :return: nothing
    """
    
    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Plotting 3D tSNE by size{bcolors.ENDC}\n")

    plt.clf()
    ax = plt.gca(projection='3d')

    fig = plt.figure(figsize=(40, 40))
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(tsne[:,0], tsne[:,1], tsne[:,2], s=50, alpha=0.6, edgecolors='w')

    ax.set_xlabel('tSNE-1')
    ax.set_ylabel('tSNE-2')
    ax.set_zlabel('tSNE-3')

    plt.tight_layout()
    plt.savefig(outpng + ".3d.png")




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='fasta input file', required=True)
    parser.add_argument('-k', help='kmer size (default=2)', default=2, type=int)
    parser.add_argument('-p', help='png output file', required=True)
    parser.add_argument('-t', help='tSNE output file')
    parser.add_argument('-d', help="dataframe tsv output file")
    parser.add_argument('-n', help='functions file. If not provided will try and extract from fasta file')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()


    if args.n:
        counts = count_pairwise_no_fn(args.f, args.k, args.v)
        fns = read_fns(args.t, args.v)
    else:
        counts, fns = count_pairwise(args.f, args.k, args.v)

    df = counts_to_df(counts, args.k, args.v)
    # Since we only make the df with non zero entries, don't need to do this
    # but should perhaps filter for correlated values
    # df = filternz(dfz, True)
    # if args.v:
    #     sys.stderr.write(f"{bcolors.BLUE}Shape before removing zeros: {dfz.shape} and after: {df.shape}\n")


    if args.d:
        with open(args.d, 'w') as out:
            df.to_csv(out, "\t")

    numcol = unique_functions(fns, args.v)
    palette = sns.color_palette("bright", numcol)
    fnar = fns_to_list(fns, df, args.v)

    tsne = create_tSNE(df, fnar, args.t, args.v)
    
    outpng = args.p.replace('.png', '')
    plot_2d(tsne, fnar, palette, outpng, args.v)

    # plot_2d_sz(tsne, fnar, palette, outpng, args.v)

    # plot_3d(tsne, fnar, palette, outpng, args.v)



