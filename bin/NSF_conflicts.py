"""
Export authors from a bib file.

You can use this to create a new conflicts list for the NSF!

Start at Google Scholar, and download all your citations as a bibtex file.

"""

from pybtex.database import parse_file
import os
import sys
import argparse
import datetime

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    now = datetime.datetime.now()
    earlyyear = now.year - 4
    parser = argparse.ArgumentParser(description='Parse a bibtex file and create a list of conflicts')
    parser.add_argument('-f', help='bibtex file', required=True)
    parser.add_argument('-c', help="Exisiting known conflicts (name, location, type)")
    parser.add_argument('-y', help="Earliest year to report conflict (default={})".format(earlyyear), default=earlyyear)
    args = parser.parse_args()

    entries = set()
    dupentries=False
    with open(args.f, 'r') as bin:
        for l in bin:
            if l.startswith('@'):
                l = l.replace('@misc', '')
                l = l.replace('@article', '')
                l = l.replace('@inproceedings', '')
                if l in entries:
                    sys.stderr.write("Duplicate entry " + l)
                    dupentries=True
                entries.add(l)

    if dupentries:
        sys.stderr.write("FATAL: The bibtex file has duplicate entries in it. Please remove them before trying to continue\n")
        sys.stderr.write("(It is an issue with Google Scholar, but pybtex breaks with duplicate entries. Sorry)\n")
        sys.exit(-1)

    bib = parse_file(args.f, 'bibtex')

    authors = set()
    authoryear = {}
    for e in bib.entries:
        if 'year' in bib.entries[e].fields:
            if int(bib.entries[e].fields['year']) > args.y:
                for aut in bib.entries[e].fields['author'].split(" and "):
                    authors.add(aut)
                    if int(bib.entries[e].fields['year']) > authoryear.get(aut, 0):
                        authoryear[aut] = int(bib.entries[e].fields['year'])

    known = {}
    if args.c:
        with open(args.c, 'r') as cin:
            for l in cin:
                l = l.strip()
                p = l.split("\t")
                known[p[1]] = l

    for a in sorted(authors):
        if a in known:
            print(known[a])
        else:
            print("C:\t{}\t{}".format(a, authoryear[a]))
