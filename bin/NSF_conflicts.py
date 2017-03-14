"""
Export authors from a bib file.

You can use this to create a new conflicts list for the NSF!

Start at Google Scholar, and download all your citations as a bibtex file.

"""

from pybtex.database import parse_file
import os
import sys
import argparse

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse a bibtex file and create a list of conflicts')
    parser.add_argument('-f', help='bibtex file', required=True)
    parser.add_argument('-c', help="Exisiting known conflicts (name, location, type)")
    args = parser.parse_args()

    entries = set()
    with open(args.f, 'r') as bin:
        for l in bin:
            if l.startswith('@'):
                l = l.replace('@misc', '')
                l = l.replace('@article', '')
                l = l.replace('@inproceedings', '')
                if l in entries:
                    sys.stderr.write("Duplicate entry " + l)
                entries.add(l)

    bib = parse_file(args.f, 'bibtex')

    authors = set()
    for e in bib.entries:
        if 'year' in bib.entries[e].fields:
            if int(bib.entries[e].fields['year']) > 2013:
                for aut in bib.entries[e].fields['author'].split(" and "):
                    authors.add(aut)

    known = {}
    if args.c:
        with open(args.c, 'r') as cin:
            for l in cin:
                l = l.strip()
                p = l.split("\t")
                if (p[3] == "Advisee" or p[3] == "Advisor"):
                    print(l)
                else:
                    known[p[1]] = l
    for a in authors:
        if a in known:
            print(known[a])
        else:
            print("Edwards, Robert\t" + a)
