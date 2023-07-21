"""
Add the function to the mmseqs easy-taxonomy file.

We want a gzip compressed version of the tophit_report file, and then we look for uniprot IDs in the first column

We use our sqlite database that has the seed to sprot and trembl connections, and the seed subsystems data
"""
import errno
import gzip
import os
import re
import sys
import argparse
from roblib import is_gzip
import sqlite3

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='mmseqs easy_taxonomy tophit_report.gz', required=True)
    parser.add_argument('-d', help='sqlite3 database of trembl, sprot, seed', default="dummy_database.sqlite")
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not os.path.exists(args.d):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.d)

    try:
        con = sqlite3.connect(args.d)
    except sqlite3.Error as e:
        sys.stderr.write(f"ERROR Connecting to database: {args.d}\n")
        sys.stderr.write(e)
        sys.exit(-1)

    urs = re.compile(r'^UniRef\d+_(\w+)')
    cur = con.cursor()

    with gzip.open(args.f, 'rt') if is_gzip(args.f) else open(args.f, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            m = urs.match(p[0])
            if m and m.group(1):
                try:
                    # cur.execute("select distinct superclass, class, subclass, subsystem_name, func from subsystems
                    # where func in (select func from trembl where uniprot = ?);", [m.group(1)])
                    cur.execute("select func from trembl where uniprot = ?", [m.group(1)])
                except sqlite3.OperationalError as e:
                    sys.stderr.write("{}".format(e))
                    sys.stderr.write("\nWhile insert on: {}\n".format(p))
                    sys.exit()

                funcs = cur.fetchone()
                if funcs:
                    func = funcs[0]
                    try:
                        cur.execute("select distinct superclass, class, subclass, subsystem_name, func from "
                                    "subsystems where func = ?", [func])
                    except sqlite3.OperationalError as e:
                        sys.stderr.write("{}".format(e))
                        sys.stderr.write("\nWhile insert on: {}\n".format(p))
                        sys.exit()
                    counts = 0
                    results = []
                    for s in (cur.fetchall()):
                        # print(">>", end="", file=sys.stderr)
                        # print("|\t|".join(s), end="", file=sys.stderr);
                        # print("<<", file=sys.stderr)
                        results.append("\t".join(list(p) + [func] + list(s)))
                        counts += 1
                    if results:
                        for r in results:
                            print(f"{r}\t{int(p[1])/counts}")
                    else:
                        # print(f"Can't find a class for {m.group(1)}", file=sys.stderr)
                        print("\t".join(p + [func, "", "", "", "", "", p[1]]))
                else:
                    print("\t".join(p + ["", "", "", "", "", "", p[1]]))

            else:
                print(f"Can't parse ID from {p[0]}", file=sys.stderr)
                print("\t".join(p + ["", "", "", "", "", "", p[1]]))

    con.close()
