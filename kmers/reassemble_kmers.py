"""
This is really a prototype for something I want to write in C. Take a list of random k-mers that overlap by 1
and reassemble them using a graph.
"""

import os
import sys
import argparse
from random import shuffle

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

def reassemble(kmers, verbose=False):
    """
    Reassemble the kmers in the array kmers into a string. the kmers must overlap by a single base
    :param kmers: an array of kmers
    :param verbose: more output
    :return: the reassembled string
    """


    # a boolean array of if we have seen a given element
    seen = [False for x in kmers]
    posn = 0
    seen[posn] = True
    string = kmers[posn]


    matched = True
    iteration = 0
    while matched:
        iteration += 1
        if verbose:
            sys.stderr.write(f"Iteration: {iteration} String: {string}\n")
        matched = False
        for i, j in enumerate(kmers):
            if seen[i]:
                continue
            # check the right end
            if string.startswith(j[1:]):
                seen[i] = True
                matched = True
                string = j[0] + string
            if string.endswith(j[:-1]):
                matched = True
                seen[i] = True
                string = string + j[-1]

    return string


if __name__ == '__main__':
    oristring = "RobIsCoolIsAGreatString"
    kmers = []
    kmerlen = 6
    for i in range(len(oristring) - kmerlen+1):
        kmers.append(oristring[i:i+kmerlen])
    shuffle(kmers)
    newstring = reassemble(kmers, False)
    print(f"Ori: {oristring}\nKmers: {kmers}\nNew: {newstring}")

    kmers = ['AGCTGC','GAGCTG','CTGCTA', 'AGAGCT', 'TAGAGC', 'GCTGCT']
    newstring2 = reassemble(kmers, False)
    print(f"Kmers: {kmers}\nNew: {newstring2}")
