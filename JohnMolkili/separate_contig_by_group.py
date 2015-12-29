import argparse
import os
import sys

import roblib

__author__ = 'Rob Edwards'

'''
Read the data file and separate the contigs into different groups.

1  - 18 -> plasma
18 - 36 -> buffy coat
37 - 54 -> csf
55 - 63 -> plasma from control

We are looking for means +/- 2std from any group that is different from all others and/or higher in non control from
plasma from control

'''


parser = argparse.ArgumentParser(description='separate out the contigs based on the group of data they are in')
parser.add_argument('-d', help='data file (all reads mapped to contigs')
args = parser.parse_args()

data = {}
headers = []
cols = {'plasma': [], 'buffy': [], 'csf': [], 'control': []}
with open(args.d, 'r') as f:
    for l in f:
        p = l.strip().split("\t")
        if not headers:
            headers = p
            for i in range(1, len(headers)):
                if int(headers[i]) <= 18:
                    cols['plasma'].append(i)
                elif int(headers[i]) <= 36:
                    cols['buffy'].append(i)
                elif int(headers[i]) <= 54:
                    cols['csf'].append(i)
                elif int(headers[i]) <= 63:
                    cols['control'].append(i)
                else:
                    sys.stderr.write("Don't understand header {} at column {}\n".format(headers[i], i))
        else:
            data[p[0]] = map(int, p[1:])

allcontigs = data.keys()
allcontigs.sort()

# calculate the mean and stdev for each group and each contig
means = {}
std = {}

for contig in allcontigs:
    means[contig] = {}
    std[contig] = {}
    for sample in cols:
        testdata = [data[contig][i] for i in cols[sample]]
        means[contig][sample] = roblib.mean(testdata)
        std[contig][sample] = 2 * roblib.stdev(testdata)

    # test the NS vs Control
    if means[contig]['plasma'] - std[contig]['plasma'] > means[contig]['control'] + std[contig]['control'] and \
        means[contig]['buffy'] - std[contig]['buffy'] > means[contig]['control'] + std[contig]['control'] and \
        means[contig]['csf'] - std[contig]['csf'] > means[contig]['control'] + std[contig]['control']:
        print("\t".join(map(str, ["ALL", contig, means[contig]['plasma'], means[contig]['buffy'], means[contig]['csf'],
                                                      means[contig]['control']])))


