
import os
import sys
import argparse
from random import random

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Subsample a file and count the frequency of different percentages')
    parser.add_argument('-f', help='percent pairwise ID file to subsample', required=True)
    parser.add_argument('-n', help="number of pairs to subsample", required=True)
    parser.add_argument('-r', help='number of times to repeat the subsample', required=True)
    parser.add_argument('-d', help='directory to write the repeat samples out to', required=True)

    args = parser.parse_args()

    if not os.path.exists(args.d):
        os.mkdir(args.d)


    lines_in_file = 100000 # for the first pass, we just use a low ball estimate of how many lines. Subsequently we'll use the real number


    for repeat in range(int(args.r)):
        sample_every  = 1.0 * int(args.n) / lines_in_file
        lines_in_file = 0
        sl = [] # sorted list
        with open(args.f, 'r') as f:
            for l in f:
                lines_in_file+=1
                if random() < sample_every:
                    p=l.strip().split("\t")
                    val = float(p[2])
                    added = False
                    for i in range(len(sl)):
                        if sl[i] > val:
                            sl[i:i]=[val]
                            added=True
                            break
                    if not added:
                        sl.append(val)
                if len(sl) >= int(args.n):
                    break

        sys.stderr.write("Repeat {} sampled the file {} times\n".format(repeat, len(sl)))

        indices = [-1 for i in range(100)]
        for i in range(len(sl)):
            ni = int(sl[i])+1
            sys.stderr.write("Adding {} for {}\n".format(ni, sl[i]))
            indices[ni]=i

        with open(os.path.join(args.d, "repeat_{}.tsv".format(repeat)), 'w') as out:
            lastindex=-1
            for i in range(50, 70):
                if indices[i] == -1:
                    indices[i] = lastindex
                else:
                    lastindex = indices[i]
                out.write("{}\t{}\n".format(i, indices[i]+1))

