"""
Generate a list of overlapping orfs from the locations.txt file
"""

import os
import sys
import argparse


locations = {}
lengths = {}
with open('locations.txt', 'r') as f:
    for l in f:
        p = l.strip().split("\t")
        if p[2] not in locations:
            locations[p[2]] = {"THEA" : [], "NONE": [], "ANY" : []}
        locations[p[2]][p[0]].append(p)
        lengths[p[2]] = int(p[6])

for contig in locations:
    posn = [[] for x in range(lengths[contig]+1)]
    # first map the NONE orfs
    for o in locations[contig]["NONE"]:
        if o[3] > o[4]:
            (o[4], o[3]) = (o[3], o[4])
        for z in range(o[3], o[4]+1):
            posn[z].append(o[1])

    # now we can find the overlaps
    seen = set()
    for o in locations[contig]["THEA"]:
        overlaps = set()
        if o[3] > o[4]:
            (o[4], o[3]) = (o[3], o[4])
        for z in range(o[3], o[4] + 1):
            overlaps.update(posn[z])
        print("\t".join(map(str, ["THEA", o[1], list(overlaps)])))
        seen.update(overlaps)

    for o in locations[contig]["ANY"]:
        overlaps = set()
        if o[3] > o[4]:
            (o[4], o[3]) = (o[3], o[4])
        for z in range(o[3], o[4] + 1):
            overlaps.update(posn[z])
        print("\t".join(map(str, ["ANY", o[1]] + list(overlaps))))
        seen.update(overlaps)

    # finally print out the NONE with no overlaps
    for o in locations[contig]["NONE"]:
        if o[1] not in seen:
            print("\t".join(map(str, ["NONE", o[1]])))

