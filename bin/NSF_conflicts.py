"""
You can use this form to update your NSF conflicts form. You'll need to do some editing, but
this will help.

First, go to Google Scholar, and visit your profile.

Select all your references by clicking the square box next to Title

Choose Export --> BibTex and then "All my articles"

That will possibly open in a new window, if so right click and choose "Save As ... " to download the file and call
it `citations.txt`

If you have a previous conflicts file, take just your conflicts and make a separate file that has five columns as
in the NSF format:
[type of conflict, conflict, university, department/email, date].

Don't worry if you dont have date, you can leave that blank. Save the file as "Tab separated text" and call that
`current_conflicts.txt`

Run this code with python3:

`python3 NSF_conflicts.py -f citations.txt -c current_conflicts.csv > revised_conflicts.tsv`

The first time you run it, it will likely complain of duplicate conflicts. That is an
issue with Google Scholar, but it breaks everything. Just edit the `citations.txt` file and remove the first
instance of each of the conflicts.

This will give you a new file called revised_conflicts.tsv that has all your conflicts. You need to go and
edit that file because you will have duplicates because some of your citations have first names and some just
have initials.


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
                    sys.stderr.write("Duplicate entry " + l.replace('{', '').replace(',', ''))
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
                if len(p) > 4:
                    if '/' in p[4]:
                        pyr = int(p[4].split('/')[2])
                    else:
                        pyr = int(p[4])

                    if p[1] in authoryear and authoryear[p[1]] and authoryear[p[1]] > pyr:
                        p[4] = "1/1/{}".format(authoryear[p[1]])
                    elif pyr and pyr > authoryear.get(p[1], 0):
                        authoryear[p[1]] = pyr
                elif p[1] in authoryear:
                    p.append("1/1/{}".format(authoryear[p[1]]))
                known[p[1]] = "\t".join(map(str, p))
                authors.add(p[1])

    for a in sorted(authors):
        toprint="" # this is just so we can fix all the unicode characters in one go
        if a in known:
            toprint = known[a]
        else:
            toprint = "C:\t{}\tUnknown\t\t1/1/{}".format(a, authoryear[a])

        toprint = toprint.replace(r'{\'e}', u"\u00E9")
        toprint = toprint.replace(r"{\~a}", u"\u00E3")
        toprint = toprint.replace(r'{\"u}', u"\u00FC")
        toprint = toprint.replace(r'{\'a}', u"\u00E1")
        toprint = toprint.replace(r"{\'o}", u"\u00F3")
        print(toprint)
