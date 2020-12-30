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

def testp(fastaf, exactp, degeneratep, verbose=False):
    for seqid, seq in stream_fasta(fastaf):
        seq = seq.upper() # always work in uppercase
        rcseq = rc(seq)
        for pp in exactp:
            lf = rf = lr = rr = 'Not Found'
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
                lr = pos+1
            except ValueError:
                pass
            try:
                pos = rcseq.index(pp[1])
                rr = pos+1
            except ValueError:
                pass
            pp += [lf, rf, lr, rr]
            print(seqid + "\t" + "\t".join(map(str, pp)))

        for pp in degeneratep:
            lre = re.compile(pp[0])
            rre = re.compile(pp[1])
            lf = rf = lr = rr = 'Not Found'
            m = lre.search(seq)
            if m:
                lf = m.pos + 1
            m = lre.search(rcseq)
            if m:
                lr = m.pos + 1
            m = rre.search(seq)
            if m:
                rf = m.pos + 1
            m = rre.search(rcseq)
            if m:
                rr = m.pos + 1
            pp += [lf, rf, lr, rr]
            print(seqid + "\t" + "\t".join(map(str, pp)))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Match primers")
    parser.add_argument('-p', help='primers file. This should be [fwd\trev]', required=True)
    parser.add_argument('-f', help='fasta filename to search', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    exactp, degeneratep =  load_primers(args.p, args.v)
    testp(args.f, exactp, degeneratep, args.v)
