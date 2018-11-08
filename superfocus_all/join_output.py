"""
Join all the output files. This may take a lot of memory :)
"""

import os
import sys
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Join parsed SF output files")
    parser.add_argument('-d', help='directory with all the output files', required=True)
    parser.add_argument('-o', help='Output file to write', required=True)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    data = {}
    sss = set()
    for infile in os.listdir(args.d):
        with open(os.path.join(args.d, infile), 'r') as f:
            l = f.readline().strip()
            if not l:
                continue

            try:
                (lvl, sample) = l.split("\t")
            except:
                sys.stderr.write("Can't split: |{}|\n".format(l))
                sys.exit(-1)
            sample = sample.replace('superfocus.bins/', '')
            data[sample] = {}
            for l in f:
                p=l.strip().split("\t")
                data[sample][p[0]] = p[1]
                sss.add(p[0])

    allss = sorted(sss)
    samples = sorted(data.keys())
    with open(args.o, 'w') as out:
        out.write("\t".join(samples))
        out.write("\n")
        for ss in allss:
            out.write(ss)
            for sam in sorted(data.keys()):
                out.write("\t{}".format(data[sam].get(ss, 0)))
            out.write("\n")
