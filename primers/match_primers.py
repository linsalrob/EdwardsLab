"""
Match primers from a list of primers to a fasta file.

Allows degenerate primers and looks in both orientations
"""

import os
import sys
import argparse
from roblib import rc, stream_fasta, colours
import re

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

def is_degenerate(seq, verbose=False):
    rep = re.compile('[AGCT]')
    res = rep.sub('', seq)
    if len(res) > 0:
        return True
    return False

def replace_degenerate(seq, verbose=False):
    oseq = seq
    seq = seq.replace('Y', '[CT]')
    seq = seq.replace('R', '[AG]')
    seq = seq.replace('S', '[GC]')
    seq = seq.replace('W', '[AT]')
    seq = seq.replace('K', '[GT]')
    seq = seq.replace('M', '[AC]')
    seq = seq.replace('B', '[GCT]')
    seq = seq.replace('D', '[AGT]')
    seq = seq.replace('H', '[ACT]')
    seq = seq.replace('V', '[ACG]')
    seq = seq.replace('N', '[ACGT]')
    rep = re.compile('[\[\]AGCT]')
    res = rep.sub('', seq)
    if len(res) > 0:
        sys.stderr.write(f"{colours.WARNING}ERROR: We replaced IUPAC characters but appear to have additional characters: {oseq} --> {seq}\n")
    return seq

def load_extended_primers(primerf, verbose=False):
    """
    Load the primers. We expect two columns, forwards and reverse
    We test for degenerate primers and if so add the regular expressions to the list
    :param primerf: the file of primers

    """
    exact = []
    degenerate = []
    names = {}

    with open(primerf, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            if (len(p) != 3):
                sys.stderr.write(f"{colours.RED}ERROR:{colours.ENDC} Can't find the primers in {l}\n")
                sys.exit(2)
            p[1] = p[1].upper()
            p[2] = p[2].upper()
            if is_degenerate(p[1]) or is_degenerate(p[2]):
                dp1 = replace_degenerate(p[1])
                degenerate.append([dp1, replace_degenerate(p[2])])
                names[dp1] = p[0]
            else:
                exact.append([p[1],p[2]])
                names[p[1]]=p[0]

    return exact, degenerate, names

def load_primers(primerf, verbose=False):
    """
    Load the primers. We expect two columns, forwards and reverse
    We test for degenerate primers and if so add the regular expressions to the list
    :param primerf: the file of primers

    """
    exact = []
    degenerate = []

    with open(primerf, 'r') as f:
        for l in f:
            l = l.upper() # always work in uppercase!
            p = l.strip().split("\t")
            if (len(p) != 2):
                sys.stderr.write(f"{colours.RED}ERROR:{colours.ENDC} Can't find the primers in {l}\n")
                sys.exit(2)
            if is_degenerate(p[0]) or is_degenerate(p[1]):
                degenerate.append([replace_degenerate(p[0]), replace_degenerate(p[1])])
            else:
                exact.append(p)

    return exact, degenerate

def testp(fastaf, exactp, degeneratep, names, printname=False, verbose=False):
    for seqid, seq in stream_fasta(fastaf):
        seq = seq.upper() # always work in uppercase
        rcseq = rc(seq)
        for pp in exactp:
            lf = rf = lr = rr = '-'
            try:
                pos = seq.index(pp[0])
                lf = pos+1
            except ValueError:
                pass
            try:
                pos = seq.index(pp[1])
                rf = pos+1
            except ValueError:
                pass
            try:
                pos = rcseq.index(pp[0])
                lr = (len(rcseq) - pos) + 1 + len(pp[0])
            except ValueError:
                pass
            try:
                pos = rcseq.index(pp[1])
                rr = (len(rcseq) - pos) + 1 + len(pp[1])
            except ValueError:
                pass
            if names:
                ttpp = [seqid, names[pp[0]]] + pp + [lf, rf, lr, rr]
            else:
                ttpp = [seqid] + pp + [lf, rf, lr, rr]

            if printname:
                print(fastaf + "\t", end="")
            print("\t".join(map(str, ttpp)))

        for pp in degeneratep:
            lre = re.compile(pp[0])
            rre = re.compile(pp[1])
            lf = rf = lr = rr = '-'
            m = lre.search(seq)
            if m:
                lf = m.pos + 1
            m = lre.search(rcseq)
            if m:
                lr = (len(rcseq) - m.pos) + 1 + len(pp[0])
            m = rre.search(seq)
            if m:
                rf = m.pos + 1
            m = rre.search(rcseq)
            if m:
                rr = (len(rcseq) - m.pos) + 1 + len(pp[1])
            if names:
                ttpp = [seqid, names[pp[0]]] + pp + [lf, rf, lr, rr]
            else:
                ttpp = [seqid] + pp + [lf, rf, lr, rr]

            if printname:
                print(fastaf + "\t", end="")
            print("\t".join(map(str, ttpp)))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Match primers")
    parser.add_argument('-x', help='extended primers file. This should be [name\tfwd\trev]')
    parser.add_argument('-p', help='primers file. This should be [fwd\trev]')
    parser.add_argument('-f', help='fasta filename to search', required=True)
    parser.add_argument('-n', help='Include file name in output', action='store_true')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()


    if args.x:
        exactp, degeneratep, names = load_extended_primers(args.x, args.v)
    elif args.p:
        exactp, degeneratep =  load_primers(args.p, args.v)
        names = None
    else:
        sys.stderr.write(f"{colours.RED}FATAL: {colours.ENDC} Either -p or -x must be supplied\n")
        sys.exit(2)


    if not os.path.exists(args.f):
        sys.stderr.write(f"{colours.RED}FATAL: {colours.ENDC} {args.f} does not exist\n")
        sys.exit(2)
    if args.v:
        print(f"Primers are\n{exactp}\n{degeneratep}")
    testp(args.f, exactp, degeneratep, names, args.n, args.v)
