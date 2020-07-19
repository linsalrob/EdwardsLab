"""
Compare bibtex output from Google and ORC ID and generate a list of papers that are in Google that are not in ORC ID
"""

import os
import sys
import argparse
from pybtex.database import parse_file, BibliographyData
from roblib import message

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


def check_dups(bibtexf, verbose=False):
    """
    Check for all duplicate entries at once.
    :param bibtexf: the bibtex file
    :param verbose: more output
    :return:
    """

    if verbose:
        message(f"Checking for duplicate entries: {bibtexf}", "PINK")
    entries = set()
    dupentries=False
    with open(bibtexf, 'r') as bin:
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


def parse_refs(bibtexf, verbose=False):
    """
    Parse the references and return some data structure
    :param bibtexf: the bibtex file
    :param verbose: more output
    :return: the BibliographyData object and a dictionary linking lower case titles to entry keys
    """

    if verbose:
        message(f"Parsing {bibtexf}", "GREEN")
    bib = parse_file(bibtexf, 'bibtex')
    titles = {}
    for e in bib.entries:
        try:
            if 'title' in bib.entries[e].fields:
                # sys.stderr.write(f"{bcolors.BLUE}{bib.entries[e].fields['title'].lower()}{bcolors.ENDC}\n")
                t = bib.entries[e].fields['title'].lower()
                t = t.replace('{', '')
                t = t.replace('}', '')
                titles[t.lower()] = e
        except Exception as ex:
            sys.stderr.write(f"Error parsing entry: {e}\n")
            print(ex)

    if verbose:
        message(f"Found {len(titles)} references", "BLUE")
    return bib, titles


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-g', help='google bibtex file', required=True)
    parser.add_argument('-o', help='orcid bibtex file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    check_dups(args.g, args.v)
    check_dups(args.o, args.v)

    g_refs, g_titles = parse_refs(args.g, args.v)
    o_refs, o_titles = parse_refs(args.o, args.v)

    print(f"Found {len(g_titles)} references in the Google bibtex and {len(o_titles)} references in the ORCID bibtex")
    dg = set(g_titles.keys()).difference(set(o_titles.keys())) # present in google but not in orcid
    do = set(o_titles.keys()).difference(set(g_titles.keys())) # present in orcid but not google
    print(f"{len(dg)} references are unique to Google\n{len(do)} references are unique to ORCID")

    g_bib = BibliographyData(entries={g_titles[e]: g_refs.entries[g_titles[e]] for e in dg})
    o_bib = BibliographyData(entries={o_titles[e]: o_refs.entries[o_titles[e]] for e in do})

    g_bib.to_file("refs_google_not_orcid.bib")
    o_bib.to_file("refs_orcid_not_google.bib")

