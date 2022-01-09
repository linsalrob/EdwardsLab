"""
Compare my titles from google sheets (abstracted to a list of just titles) to the references in paperpile
available

NOTE: SEE https://github.com/linsalrob/CompareReferences
"""

import os
import sys
import argparse

from roblib import bcolors
from roblib import bcolors
from pybtex.database import parse_file

doctitles = set()
with open("/home/redwards/Desktop/MyRefs/my_titles.txt", "r") as f:
    for l in f:
        l = l.strip()
        if l.endswith("."):
            l = l.rsplit(".", 1)[0]
        if l in doctitles:
            sys.stderr.write(f"{bcolors.RED}Already have {l}{bcolors.ENDC}\n")
            continue
        doctitles.add(l.lower())

bib = parse_file("/home/redwards/Desktop/MyRefs/paperpile.bib", 'bibtex')

bibtitles = set()
for e in bib.entries:
    try:
        if 'title' in bib.entries[e].fields:
            # sys.stderr.write(f"{bcolors.BLUE}{bib.entries[e].fields['title'].lower()}{bcolors.ENDC}\n")
            t = bib.entries[e].fields['title'].lower()
            t = t.replace('{', '')
            t = t.replace('}', '')
            bibtitles.add(t)
    except Exception as ex:
        sys.stderr.write(f"Error parsing entry: {e}\n")
        print(ex)


# index the first n characters
n = 12
docindex = {}
doctindex = {}
for t in doctitles:
    docindex[t[:n]] = t
    doctindex[t] = t[:n]
bibindex = {}
bibtindex = {}
for t in bibtitles:
    bibindex[t[:n]] = t
    bibtindex[t] = t[:n]

if False:
    bibnotdoc = bibtitles.difference(doctitles)
    # many of these differences are not real.
    alttitles = []
    for t in bibnotdoc:
        if bibtindex[t] in docindex:
            alttitles.append((t, docindex[bibtindex[t]]))
        else:
            print(t)
else:
    bibnotdoc = doctitles.difference(bibtitles)
    # many of these differences are not real.
    alttitles = []
    for t in bibnotdoc:
        if doctindex[t] in bibindex:
            alttitles.append((t, bibindex[doctindex[t]]))
        else:
            print(t)

print("\n\nALTERNATE TITLES\n")
for a in alttitles:
    print(f"{a[0]}\n{a[1]}\n")

exit(0)
for t in doctitles:
    if t.startswith('microfluidic'):
        print(t)
print("\n")
for t in bibtitles:
    if t.startswith('microfluidic'):
        print(t)

