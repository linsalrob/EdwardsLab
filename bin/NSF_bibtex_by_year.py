"""
Parse a bibtex file and only print those entries within a certain number of years

"""

import os
import sys
import argparse
import datetime
from pybtex.database import parse_file, BibliographyData

if __name__ == "__main__":
    now = datetime.datetime.now()
    earlyyear = now.year - 4
    parser = argparse.ArgumentParser(description='Parse a bibtex file and create a list of conflicts')
    parser.add_argument('-f', help='bibtex file', required=True)
    parser.add_argument('-y', help="Earliest year to report conflict (default={})".format(earlyyear), default=earlyyear, type=int)
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
                    sys.stderr.write("Duplicate entry " + l.replace('{', '').replace(',', ''))
                    dupentries=True
                entries.add(l)

    if dupentries:
        sys.stderr.write("FATAL: The bibtex file has duplicate entries in it. Please remove them before trying to continue\n")
        sys.stderr.write("(It is an issue with Google Scholar, but pybtex breaks with duplicate entries. Sorry)\n")
        sys.exit(-1)

    bib = parse_file(args.f, 'bibtex')

    for e in bib.entries:
        if 'year' in bib.entries[e].fields:
            if int(bib.entries[e].fields['year']) >= args.y:
                bib_data = BibliographyData({e : bib.entries[e]})
                print(bib_data.to_string('bibtex'))


