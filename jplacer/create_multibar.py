"""
Create a multibar file from our labeled file and our tree. We will do this at different taxonomic levels
and different label levels.

For each label level we have, we'll create a single file.

"""

import os
import sys
import argparse
from ete3 import Tree

def read_labels(lf, col, verbose=False):
    """
    Read the labels file and return a dict with tree labels and values
    :param lf: labels file
    :param col: the column to use
    :param verbose: extra output
    :return: a dict of the leaves and their labels and a dict of the labels and their counts
    """

    ret = {}
    counts = {}
    with open(lf, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            if len(p) < col:
                continue
            try:
                if not p[col]:
                    continue
            except:
                sys.stderr.write("Error with: {} and col: {}\n".format(l.strip(), col))
                continue
            ret[p[0]] = p[col]
            counts[p[col]] = counts.get(p[col], 0) + 1
    return ret, counts



def write_directory(treefile, data, counts, taxa, outputdir, colors, legend, legendshape, proportions, verbose=False):
    """
    Write a directory with one multibar file per type
    :param treefile: The tree file to parse
    :param data: The data dict with leaves and labels
    :param counts: the counts of label frequency
    :param outputdir: the directory to create
    :param colors: the array of colors to choose from
    :param legend: the legend name
    :param legendshape: the legend shape
    :param verbose: more output
    :return:
    """


    allkeys = list(counts.keys())
    if len(allkeys) > len(colors):
        sys.stderr.write("ERROR: Not enough colors. We have {}  keys and {} colors\n".format(len(allkeys), len(colors)))
        sys.exit(-1)
    keycolors = {x: colors[allkeys.index(x)] for x in allkeys}

    if not os.path.exists(outputdir):
        try:
            os.mkdir(outputdir)
        except Exception as e:
            sys.stderr.write("Cannot make directory: {}\n".format(outputdir))
            sys.stderr.write("{}\n".format(e))
            sys.exit(-1)
    if verbose:
        sys.stderr.write("Reading tree\n")

    tree = Tree(treefile, quoted_node_names=True, format=1)

    if verbose:
        sys.stderr.write(f"Creating output files in {outputdir}\n")

    for k in counts:
        outputf = os.path.join(outputdir, k + ".multibar.txt")

        with open(outputf, 'w') as out:
            out.write("DATASET_MULTIBAR\nSEPARATOR COMMA\n")
            out.write("DATASET_LABEL,{} counts\n".format(k))
            out.write("FIELD_COLORS,{}\n".format(keycolors[k]))
            out.write("FIELD_LABELS,{}\n".format(k))
            out.write("WIDTH,50\n")
            out.write("DATASET_SCALE,0-{}-{}\n".format(k, keycolors[k]))
            out.write("HEIGHT_FACTOR,50\n")
            out.write("SHOW_INTERNAL,1\n")
            out.write("ALIGN_FIELDS,1\n")
            out.write("COLOR,{}\n".format(keycolors[k]))
            out.write("DATA\n")

            for n in tree.traverse("preorder"):
                if taxa in n.name:
                    leafcount = 0
                    for l in n.get_leaves():
                        if l.name in data and data[l.name] == k:
                            leafcount += 1
                    if proportions:
                        leafcount /= counts[k]
                    out.write("{},{}\n".format(n.name, leafcount))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-f', help='The labeled leaves file from fastq2ids.py', required=True)
    parser.add_argument('-t', help='Newick tree file', required=True)
    parser.add_argument('-d', help='Output directory where to write the files', required=True)
    parser.add_argument('-n', help='Column in the labeled leaves file to use. 0 indexed', required=True, type=int)
    parser.add_argument('-l', help='Color strip legend (e.g. Kingdom, Fish, Species', required=True)
    parser.add_argument('-x', help='taxa to use for the labels', required=True)
    parser.add_argument('-s', help='Legend shape (a number). Default = 1', default="1", type=str)
    parser.add_argument('-p', help='Display proportion of counts not counts', action='store_true')
    parser.add_argument('-c', help='Colors to use. These will be prepended to our default list', action='append')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999']
    if args.c:
        colors = args.c + colors

    data, counts = read_labels(args.f, args.n, args.v)

    taxa= args.x
    if not taxa.startswith('r_'):
        taxa = "r_{}".format(taxa)

    write_directory(args.t, data, counts, taxa, args.d, colors, args.l, args.s, args.p, args.v)