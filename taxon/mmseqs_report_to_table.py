"""
Convert a directory of directories of mmseqs outputs to a directory of tables
"""

import os
import sys
import argparse

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-d', help='directory of directories', required=True)
    parser.add_argument('-o', help='output directory', required=True)
    parser.add_argument('-e', help='file extension. [Default %(default)s]', default='_report.gz')
    parser.add_argument('-i', help='file extension to ignore. [Default %(default)s]', default='_tophit_report.gz')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    data = {}
    wanted = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    wanted_lut = {x:wanted.index(x) for x in wanted}

    allsamples = {"all"}
    allnames = {w:set() for w in wanted}

    for d in os.listdir(args.d):
        for f in os.listdir(os.path.join(args.d, d)):
            if f.endswith(args.e) and not f.endswith(args.i):
                allsamples.add(f)
                with open(os.path.join(args.d, d, f), 'r') as rep:
                    data[f] = {"all" : {}}
                    current = ['', '', '', '', '', '', '']
                    for l in rep:
                        p = l.strip().split("\t")
                        taxa = p[3]
                        if taxa not in wanted_lut:
                            continue
                        name = f"{taxa[0]}__{p[4].strip()}"
                        allnames[taxa].add(name)
                        count = p[1]
                        if taxa not in data[f]:
                            data[f][taxa] = {}
                        data[f][taxa][name] = count
                        if taxa == 'species':
                            current[6] = name
                            data[f]["all"][";".join(current)] = count
                        else:
                            i = wanted_lut[taxa]
                            current[i] = name
                            for j in range(i+1, 7):
                                current[j] = ""


    # open output files for straight counts
    samples = sorted(list(allsamples))
    sample_names = "\t".join(samples)
    for w in wanted:
        filename = f"{w}.tsv"
        if w == 'all':
            filename = "all_levels.tsv"
        with open(os.path.join(args.o, filename), 'w') as out:
            print(f"#NAME\t{sample_names}", file=out)
            for n in allnames[w]:
                print(n, file=out, end="")
                for s in samples:
                    if n in data[s]['all']:
                        print(f"\t{data[s][w][n]}", file=out, end="")
                    else:
                        print("\t0", file=out, end="")
                print(file=out)
