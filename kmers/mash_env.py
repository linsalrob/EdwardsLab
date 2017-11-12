"""
Print out all combinations of mash distances and environments
"""

import os
import sys
import argparse
import gzip

def read_env(f):
    """
    Read the environment list. This is a file with genome id\tEnvironment
    :param f:
    :return:
    """
    env = {}
    with open(f, 'r') as i:
        for l in i:
            p=l.strip().split("\t")
            env[p[0]] = p[1]
    return env

def read_mash(f, env):
    """
    Read the mash output file. This is the massive file from bob that has all pairwise mash scores
    :param f: the mash output file
    :param env:  the environmant hash from read_env
    :return:
    """

    i=gzip.open(f)
    for l in i:
        p=l.strip().split("\t")
        g1 = p[0].split("/")[-1].replace(".fna", '')
        g2 = p[1].split("/")[-1].replace(".fna", '')
        dist = p[2]
        if g1 in env and g2 in env:
            print("{}\t{}\t{}".format(env[g1], env[g2], dist))




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-f', help='mash file', required=True)
    parser.add_argument('-e', help='environments file', required=True)
    parser.add_argument('-v', help='verbose output')
    args = parser.parse_args()

    env = read_env(args.e)
    read_mash(args.f, env)
    