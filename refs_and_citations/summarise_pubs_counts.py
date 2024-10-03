"""
Once you have (finally!!) updated the number of citations in the excel file, this script will summarise the
data based on the x's in the columns 25 onwards.
"""

import os
import sys
import argparse

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='My Publications file', required=True)
    parser.add_argument('-c', help='citations columns', type=int, required=True)
    parser.add_argument('-b', help="first column from where to count x's", default=23, type=int)
    parser.add_argument('--merge', help='merge these two (or more) columns', type=int, action='append')
    parser.add_argument('--merge_and', help='merge with an AND not an OR', action='store_true')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()


    header = []
    counts = {}
    combined_header = None
    with open(args.f, 'r') as file:
        for l in file:
            p = l.rstrip().split("\t")
            if not header:
                header = p
                for i in range(args.b,len(p)):
                    counts[header[i]] = []
                if args.merge:
                    newhead = []
                    for c in args.merge:
                        newhead.append(header[c])
                    if args.merge_and:
                        combined_header = " AND ".join(newhead)
                    else:
                        combined_header = " OR ".join(newhead)
                    counts[combined_header] = []
                continue
            for i in range(args.b, len(p)):
                if 'x' in p[i]:
                    if p[args.c]:
                        counts[header[i]].append(int(p[args.c]))
            if args.merge:
                newsum = 0
                nin = 0
                for c in args.merge:
                    if c >= len(p):
                        continue
                    if 'x' in p[c]:
                        if p[args.c]:
                            newsum += int(p[args.c])
                            nin += 1
                if newsum > 0:
                    if args.merge_and:
                        if nin == len(args.merge):
                            counts[combined_header].append(newsum)
                    else:
                        counts[combined_header].append(newsum)

    for label in counts:
        print(f"{label} : {len(counts[label])} papers; {sum(counts[label])} citations")

