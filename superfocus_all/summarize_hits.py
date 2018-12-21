"""
Summarize the hits from the superfocus all analysis. We normalize the hits to the
length of the metagenome and summarize by subsystem
"""

import os
import sys
import argparse
import gzip

def read_counts(rdir, verbose=False):
    """
    Read the counts file and return a dict of sequence and total number of bp
    :param rdir: the results directory
    :param verbose: more output
    :return: a dict of sample -> number of bp
    """

    if not os.path.exists(os.path.join(rdir, "counts.txt.gz")):
        sys.stderr.write("Error: {} does not exist. Can not read counts\n".format(os.path.join(rdir, "counts.txt.gz")))
        sys.exit(-1)

    seqlen = {}
    with gzip.open(os.path.join(rdir, "counts.txt.gz"), 'rt') as f:
        seq1 = f.readline().strip().split("/")[-1]
        n = f.readline()
        t = f.readline().split()[-1]
        seqlen[seq1] = int(t)
        while n and '--' not in n:
            n = f.readline()
        n = f.readline()
        if n:
            seq2 = n.strip().split("/")[-1]
            n = f.readline()
            t = f.readline().split()[-1]
            seqlen[seq2] = int(t)

    return seqlen


def read_bins(rdir, lvl, verbose=False):
    """
    Read and summarize the bins
    :param rdir: the results directory
    :param lvl: the level to store
    :param verbose: more output
    :return: a dict of the lvl and number of base pairs that match
    """


    binningfile = None
    for f in os.listdir(rdir):
        if f.endswith("binning.xls.gz"):
            binningfile = f

    if not binningfile:
        sys.stderr.write("ERROR: Can not find a binning file in {}\n".format(rdir))
        sys.exit(-1)

    # note that Geni had a bug in superfocus where some lines were duplicated. This eliminates that bug ;)
    results = set()
    with gzip.open(os.path.join(rdir, binningfile), 'rt') as f:
        l = f.readline()
        while l and not l.startswith('Sample name'):
            l = f.readline()
        for l in f:
            results.add(l.strip())

    if verbose:
        sys.stderr.write("Found {} unique results\n".format(len(results)))

    counts = {}
    for r in results:
        p = r.strip().split("\t")
        if p[0] not in counts:
            counts[p[0]] = {}
        thislvl = p[lvl]
        if 3 == lvl:
            # confusingly, this is level 2!
            thislvl = p[2] + " :: " + p[3]
        counts[p[0]][thislvl] = counts[p[0]].get(thislvl, 0) + float(p[7])

    if verbose:
        for s in counts:
            sys.stderr.write("Found {} unique subsystem levels in column {} for sample {}\n".format(len(counts[s]), lvl, s))

    return counts



if __name__ == '__main__':
    lvlchoices = ['1', '2', '3', 'fn']
    lvlchoicesstr = "'".join(lvlchoices)
    parser = argparse.ArgumentParser(description="Summarize the superfocus output")
    parser.add_argument('-d', help='Directory with all the results', required=True)
    parser.add_argument('-l', help='Subsystem level. Must be one of {}. Default = 1'.format(lvlchoicesstr), default='1')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    seqlen = read_counts(args.d, args.v)

    if args.l not in lvlchoices:
        sys.stderr.write("ERROR: Level must be one of '{}'\n".format("'".join(lvlchoicesstr)))
        sys.exit(-1)
    lvl = lvlchoices.index(args.l)
    lvl += 2

    counts = read_bins(args.d, lvl, args.v)

    # if we have two read files (pass_1 and pass_2), lets average them
    # if we only have one, we don't need to!

    sys.stdout.write("Subsystem {}\t{}\n".format(args.l, args.d))
    if len(counts) == 2:
        avcounts = {}
        (k1,k2) = counts.keys()
        for s in counts[k1]:
            if s not in counts[k2]:
                # sys.stderr.write("WARNING: Only found {} in one of the read files. So we use that number, and don't average it. This is probably wrong!\n".format(s))
                avcounts[s] = counts[k1][s]
            else:
                avcounts[s] = ((1.0 * counts[k1][s] / seqlen[k1]) + (1.0 * counts[k2][s] / seqlen[k2])) / 2
        for s in counts[k2]:
            if s not in counts[k1]:
                # sys.stderr.write("WARNING: Only found {} in one of the read files. So we use that number, and don't average it. This is probably wrong!\n".format(s))
                avcounts[s] = counts[k2][s]
        for s in avcounts:
            sys.stdout.write("{}\t{}\n".format(s, avcounts[s]))
    elif len(counts) == 1:
        k1 = list(counts.keys())[0]
        for s in counts[k1]:
            sys.stdout.write("{}\t{}\n".format(s, (1.0 * counts[k1][s] / seqlen[k1])))
    else:
        sys.stderr.write("Bugger. We found {} keys in counts. Don't know what to do!\n".format(len(counts)))
