"""
Count metagenomes under phylum level
"""

import os
import sys
import argparse
from ete3 import Tree


def read_mg_types(mgfile):
    """
    Read the metagenome types from a file where the first four columns
    are [domain, type, species, seqid]. For example, sharks_fish.distance.labelled.tsv
    which is the output from join.pl
    :param mgfile: the sharks_fish.distance.labelled.tsv file
    :return: a dict of all the metagenomes
    """
    
    val = {}
    with open(mgfile, 'r') as f:
        for l in f:
            p=l.strip().split("\t")
            if "Metagenome" in p[0]:
                val[p[3]] = [p[1], p[2]]
    return val

def write_one_file(tree, taxa, fish, fishcount, typecount, types, datasetcounts, c, outputf):

    with open(outputf, 'w') as out:
        out.write("DATASET_MULTIBAR\nSEPARATOR COMMA\n")
        out.write("DATASET_LABEL,Metagenome {} counts\n".format(taxa.replace("r_", "")))
        out.write("SHOW_INTERNAL,1\n")
        out.write("ALIGN_FIELDS,1\n")
        out.write("COLOR,#42F4CB\n")
        out.write("FIELD_COLORS,{}\n".format(",".join([c[x] for x in fish + types])))
        out.write("FIELD_LABELS,{}\n".format(",".join(fish + types)))
        out.write("LEGEND_TITLE,Metagenome Counts\n")
        out.write("LEGEND_SHAPES,{},{}\n".format(",".join("1" for x in fish), ",".join("2" for x in types)))
        out.write("LEGEND_COLORS,{}\n".format(",".join([c[x] for x in fish + types])))
        out.write("LEGEND_LABELS,{}\n".format(",".join(fish + types)))
        out.write("DATA\n")
        """
        Note: THis block writes the dataset for only the nodes of interest.
        Below we have a hack to write the dataset for all nodes :)

        for p in fishcount:
            out.write(p)
            for f in fish:
                out.write(",{}".format(fishcount[p].get(f, 0)))
            for t in types:
                out.write(",{}".format(typecount[p].get(t, 0)))
            out.write("\n")
        """
        
        printstring = None
        printed=set()
        for n in tree.traverse("preorder"):
            if taxa in n.name:
                printstring = ""
                for f in fish:
                    printstring += ",{}".format(fishcount[n.name].get(f, 0))
                for t in types:
                    printstring += ",{}".format(typecount[n.name].get(t, 0))
                for l in n.traverse("preorder"):
                    if l.name in printed:
                        continue
                    out.write("{}{}\n".format(l.name, printstring))
                    printed.add(l.name)
                printstring = None

    


def write_directory(tree, taxa, fish, fishcount, typecount, types, datasetcounts, c, outputdir):

    if not os.path.exists(outputdir):
        try:
            os.mkdir(outputdir)
        except Exception as e:
            sys.stderr.write("Cannot make directory: {}\n".format(outputdir))
            sys.stderr.write("{}\n".format(e))
            sys.exit(-1)

    for f in fish:
        outputf = os.path.join(outputdir, f + ".multibar.txt")

        with open(outputf, 'w') as out:
            out.write("DATASET_MULTIBAR\nSEPARATOR COMMA\n")
            out.write("DATASET_LABEL,{} counts\n".format(f))
            out.write("FIELD_COLORS,{}\n".format(c[f]))
            out.write("FIELD_LABELS,{}\n".format(f))
            out.write("WIDTH,50\n")
            out.write("DATASET_SCALE,0-{}-{}\n".format(f, c[f]))
            out.write("HEIGHT_FACTOR,50\n")
            out.write("SHOW_INTERNAL,1\n")
            out.write("ALIGN_FIELDS,1\n")
            out.write("COLOR,{}\n".format(c[f]))
            out.write("DATA\n")
            
            printstring = None
            printed=set()
            for n in tree.traverse("preorder"):
                if taxa in n.name:
                    printstring = ",{}".format(fishcount[n.name].get(f, 0))
                    for l in n.traverse("preorder"):
                        if l.name in printed:
                            continue
                        out.write("{}{}\n".format(l.name, printstring))
                        printed.add(l.name)
                    printstring = None

    
    for t in types:
        outputf = os.path.join(outputdir, t + ".multibar.txt")

        with open(outputf, 'w') as out:
            out.write("DATASET_MULTIBAR\nSEPARATOR COMMA\n")
            out.write("DATASET_LABEL,{} counts\n".format(t))
            out.write("FIELD_COLORS,{}\n".format(c[t]))
            out.write("FIELD_LABELS,{}\n".format(t))
            out.write("DATASET_SCALE,0-{}-{}\n".format(t, c[t]))
            out.write("WIDTH,50\n")
            out.write("HEIGHT_FACTOR,50\n")
            out.write("SHOW_INTERNAL,1\n")
            out.write("ALIGN_FIELDS,1\n")
            out.write("COLOR,{}\n".format(c[t]))
            out.write("DATA\n")
            
            printstring = None
            printed=set()
            for n in tree.traverse("preorder"):
                if taxa in n.name:
                    printstring = ",{}".format(typecount[n.name].get(t, 0))
                    for l in n.traverse("preorder"):
                        if l.name in printed:
                            continue
                        out.write("{}{}\n".format(l.name, printstring))
                        printed.add(l.name)
                    printstring = None

    


def count_mg(treefile, mgdesc, taxa, propc=False, outputf=None, outputdir=None):
    """
    Read the tree and count below phylum
    :param treefile: the file with the tree
    :param mgdesc: the description of the metagenomes from read_mg_types
    :param taxa: the taxonomic level (phylum, class, species, etc)
    :param outputf: the output file name
    :param outputdir: write the output to a directory
    :param propc: use proportional counts not raw counts
    """

    tree = Tree(treefile, quoted_node_names=True, format=1)

    if not taxa.startswith('r_'):
        taxa = "r_{}".format(taxa)
    
    fishcount = {}
    allfish = set()
    typecount = {}
    alltypes = set()
    datasetcounts = {}
    maxval = 0
    for n in tree.traverse("preorder"):
        if taxa in n.name:
            if n.name not in fishcount:
                fishcount[n.name]={}
            if n.name not in typecount:
                typecount[n.name]={}
            for l in n.get_leaves():
                if l.name in mgdesc:
                    vals = mgdesc[l.name]
                    fishcount[n.name][vals[0]] = fishcount[n.name].get(vals[0], 0) + 1
                    if fishcount[n.name][vals[0]] > maxval:
                        maxval = fishcount[n.name][vals[0]]
                    typecount[n.name][vals[1]] = typecount[n.name].get(vals[1], 0) + 1
                    if typecount[n.name][vals[1]] > maxval:
                        maxval = typecount[n.name][vals[1]]
                    allfish.add(vals[0])
                    alltypes.add(vals[1])
                    # this assumes that the names of the types (fish/shark) and the species of fish/shark are different
                    datasetcounts[vals[0]] = datasetcounts.get(vals[0], 0)+1
                    datasetcounts[vals[1]] = datasetcounts.get(vals[1], 0)+1
    if propc:
        for n in fishcount:
            for f in fishcount[n]:
                fishcount[n][f] = fishcount[n][f]/datasetcounts[f]
        for n in typecount:
            for t in typecount[n]:
                typecount[n][t] = typecount[n][t]/datasetcounts[t]

    colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']
    fish = list(allfish)
    types = list(alltypes)
    c = {}
    for i,j in enumerate(fish + types):
        c[j] = colors[i]


    sys.stderr.write("We have {} fish and {} types\n".format(len(fish), len(types)))

    if outputf:
        write_one_file(tree, taxa, fish, fishcount, typecount, types, datasetcounts, c, outputf)
    else:
        write_directory(tree, taxa, fish, fishcount, typecount, types, datasetcounts, c, outputdir)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Count metagenomes on the tree')
    parser.add_argument('-t', help='tree file', required=True)
    parser.add_argument('-d', help='metagenome description file', required=True)
    parser.add_argument('-l', help='taxonomic level', required=True)
    parser.add_argument('-o', help='output filename')
    parser.add_argument('-r', help='output directory (one file per type)')
    parser.add_argument('-p', help='Display proportion of counts not counts', action='store_true')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if args.o and args.d:
        sys.stderr.write("Please only specify output file OR output directory, not both\n")
        sys.exit(-1)
    if not args.o and not args.d:
        sys.stderr.write("Please specify one of output directory or output file (but not both)\n")
        sys.exit(-1)

    mt = read_mg_types(args.d)
    
    if args.o:
        count_mg(args.t, mt, args.l, args.p, outputf=args.o)
    if args.d:
        count_mg(args.t, mt, args.l, args.p, outputdir=args.r)


