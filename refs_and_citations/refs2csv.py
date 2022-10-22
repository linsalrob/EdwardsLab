"""
convert a bibtex file to tsv
"""

import os
import sys
import argparse
from roblib_tk import choose_a_file, write_a_file
from bibtex import parse_bibtex_file

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='input bibtext file')
    parser.add_argument('-o', help='output tsv file')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if args.f:
        bt = parse_bibtex_file(args.f, args.v)
    else:
        filename = choose_a_file("Choose a bibtex file")
        bt = parse_bibtex_file(filename, False)

    if args.o:
        out = open(args.o, "w", encoding="utf-8")
    else:
        filename = write_a_file("Where to save the results")
        out = open(filename, "w", encoding="utf-8")



    fields = ["year", "title", "booktitle", "journal", "volume", "number", "organization", "pages", "publisher", "school", "institution"]
    persons = ["author", "editor"]

    out.write('key\tType\tAuthors\tEditors\t')
    out.write("\t".join(fields))
    out.write("\n")

    for e in bt.entries:
        data = [bt.entries[e].key, bt.entries[e].type, "", ""]
        data += ["" for w in fields]
        for f in bt.entries[e].fields:
            if f not in fields:
                sys.stderr.write(f"Found new field: '{f}'\n")
                continue
            data[fields.index(f)+4] = bt.entries[e].fields[f]

        if hasattr(bt.entries[e], 'persons'):
            if 'author' in bt.entries[e].persons:
                authors = []
                for a in bt.entries[e].persons['author']:
                    authors.append(str(a))
                data[2] = ", ".join(authors)
            if 'editor' in bt.entries[e].persons:
                editors = []
                for e in bt.entries[e].persons['editor']:
                    editors.append(str(e))
                data[3] = ", ".join(editors)
            for k in bt.entries[e].persons:
                if k not in ['author', 'editor']:
                    sys.stderr.write(f"Have a person {k} in {bt.entries[e].key}\n")
        else:
            sys.stderr.write(f"No persons in {bt.entries[e].key}\n")
        print("\t".join(data), file=out)


    out.close()
