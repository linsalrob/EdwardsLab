



import os
import sys
import argparse
import re
from roblib import stream_fastq
from taxon import get_taxonomy_db, get_taxonomy

sys.stderr.write("Connecting to db\n")
c = get_taxonomy_db() 
sys.stderr.write("Connected\n")





def determine_phylogeny(fn, verbose=False):

    m = re.search('\[(\d+)\]', fn)
    if not m:
        return "Unknown", None

    tid = m.groups()[0]
    t,n = get_taxonomy(tid, c)
    if not t:
        if verbose:
            sys.stderr.write("Can't find tax for {} in the db\n".format(tid))
        return "Unknown", None

    while t.parent > 1 and t.parent != 131567:
        # 131567 is cellular organisms
        t,n = get_taxonomy(t.parent, c)
    return n.scientific_name, None


def get_dom(leaff, verbose=False):
    with open(leaff, 'r') as f:
        for l in f:
            l = l.strip()
            dom, nm = determine_phylogeny(l, verbose)
            print("{}\t{}\t{}".format(l, dom, nm))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Test getting domains')
    parser.add_argument('-l', help='leaves list file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    get_dom(args.l, args.v)


