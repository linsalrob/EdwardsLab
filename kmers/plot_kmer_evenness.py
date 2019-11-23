"""
Plot the shannon and evennes of the kmer profiles. Uses the output from kmer_entropy.py code, but
requires that you add headers.
"""

import os
import sys
import argparse

import os, sys
import numpy as np
import pandas as pd
import seaborn as sns
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.set_style("white")
import matplotlib.pyplot as plt
from roblib import bcolors

def check_outputf(f, verbose=True):
    """
    Check the header
    :param f:
    :param verbose:
    :return:
    """


    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Checking data{bcolors.ENDC}\n")
    with open(f, 'r') as fin:
        l = fin.readline()
        if not l.startswith('file'):
            sys.stderr.write(f"{bcolors.FATAL}ERROR. Please add the header line to {f}\n")
            sys.stderr.write("It should have [file, kmer, Shannon, Evenness] separated by tabs\n{bcolors.ENDC}\n")
            sys.exit(-1)


def read_df(f, sample=None, verbose=False):

    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Reading data{bcolors.ENDC}\n")
    dfa = pd.read_csv(f, delimiter="\t")
    o = dfa.shape
    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Removing outliers\n{bcolors.ENDC}")
    dfa = dfa[dfa['file'] != "fasta/AH004327.fasta"]
    dfa = dfa[dfa['file'] != 'fasta/AB830321.fasta']

    if sample:
        dfa = dfa.sample(n=sample)
        dftest = dfa[dfa['kmer'] == 77]
        sys.stderr.write(f"{bcolors.BLUE}Sampled {sample} entries. We have {dftest.shape[0]} genomes\n")

    n = dfa.shape
    print(f"{bcolors.GREEN}Original dataframe was {o[0]} entries. After filtering we have {n[0]} entries\n{bcolors.ENDC}")
    return dfa


def plot_evenness(df, output, verbose=False):
    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Plotting evenness{bcolors.ENDC}\n")
    sns.violinplot(data=df, x='kmer', y='Evenness')
    sns.despine(offset=10, trim=True)
    plt.savefig(f"{output}.evenness.png")
    plt.clf()

def plot_swarm_evenness(df, output, verbose=False):
    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Plotting swarmed evenness{bcolors.ENDC}\n")
    sns.violinplot(data=df, x='kmer', y='Evenness')
    sns.swarmplot(data=df, x='kmer', y='Evenness')
    sns.despine(offset=10, trim=True)
    plt.savefig(f"{output}.swarm.evenness.png")
    plt.clf()

def plot_shannon(df, output, verbose=False):
    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Plotting swarmed shannon{bcolors.ENDC}\n")
    sns.violinplot(data=df, x='kmer', y='Shannon')
    sns.swarmplot(data=df, x='kmer', y='Shannon')
    sns.despine(offset=10, trim=True)
    plt.savefig(f"{output}.shannon.png")
    plt.clf()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot a heatmap")
    parser.add_argument('-f', help='kmer_entropy.py file', required=True)
    parser.add_argument('-o', help='output image file basename', required=True)
    parser.add_argument('-s', help='subsample data to make quicker smaller plots')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    check_outputf(args.f, args.v)

    df = read_df(args.f, args.s, args.v)
    output = args.o.replace('.png', '')
    plot_evenness(df, output, args.v)
    plot_swarm_evenness(df, output, args.v)
    plot_shannon(df, output, args.v)