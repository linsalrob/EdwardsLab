"""
Convert a directory of directories of mmseqs outputs to a directory of tables
"""

import os
import sys
import argparse
import gzip

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

    allsamples = set()
    allnames = {w:set() for w in wanted}
    allnames['all'] = set()

    for d in os.listdir(args.d):
        for f in os.listdir(os.path.join(args.d, d)):
            if f.endswith(args.e) and not f.endswith(args.i):
                samplename = f.replace(args.e, "")
                allsamples.add(samplename)
                if args.v:
                    print(f"Parsing {f}", file=sys.stderr)
                with gzip.open(os.path.join(args.d, d, f), 'rt') as rep:
                    data[samplename] = {}
                    data[samplename]["all"] = {}
                    current = ['', '', '', '', '', '', '']
                    for l in rep:
                        p = l.strip().split("\t")
                        taxa = p[3]
                        if taxa not in wanted_lut:
                            continue
                        name = f"{taxa[0]}__{p[5].strip()}"
                        allnames[taxa].add(name)
                        count = p[1]
                        if taxa not in data[samplename]:
                            data[samplename][taxa] = {}
                        data[samplename][taxa][name] = count
                        if taxa == 'species':
                            current[6] = name
                            ct = ";".join(current)
                            data[samplename]["all"][ct] = count
                            allnames['all'].add(ct)
                        else:
                            i = wanted_lut[taxa]
                            current[i] = name
                            for j in range(i+1, 7):
                                current[j] = ""


    # open output files for straight counts
    samples = sorted(list(allsamples))
    sample_names = "\t".join(samples)
    os.makedirs(args.o, exist_ok=True) 
    for w in wanted + ["all"]:
        filename = f"{w}.tsv"
        if w == 'all':
            filename = "all_levels.tsv"
        if args.v:
            print(f"Writing {filename}", file=sys.stderr)
        with open(os.path.join(args.o, filename), 'w') as out:
            print(f"#NAME\t{sample_names}", file=out)
            for n in allnames[w]:
                print(n, file=out, end="")
                for s in samples:
                    if 'all' not in data[s]:
                        print(f"No 'all' in {s}", file=sys.stderr)
                        continue
                    if w in data[s] and n in data[s][w]:
                        print(f"\t{data[s][w][n]}", file=out, end="")
                    else:
                        print("\t0", file=out, end="")
                print(file=out)
