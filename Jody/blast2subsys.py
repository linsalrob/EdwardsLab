"""
Read the blast (well, diamond) output and convert it to subsystem information

This uses a sql database I made for the protein IDs and subsystems. The database
is in two locations:
- /home3/redwards/Jody/Subsystems2020/protein_functions_subsystems.sql  ## the original database
- /home/redwards/protein_functions_subsystems.sql ## on rambox only, a local copy on spinning disk so it doesn't suffer nfs issues

 sqlite3 -separator $'\t' /home3/redwards/Jody/Subsystems2020/protein_functions_subsystems.sql 'select * from proteins inner join md5sum on proteins.patric_id = md5sum.patric_id where md5sum.md5sum = "2755d47229436e24f7548ea7a2bd6443";'
"""

import os
import sys
import argparse
import sqlite3
from roblib import message, stream_blast_results

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

proteins = {}
roles = {}

def connect_to_db(dbname, verbose=False):
    """
    Connect to the database
    :param dbname: the database file name
    :param verbose: print addtional output
    :return: the database connection
    """

    global conn
    if conn:
        return conn

    try:
        conn = sqlite3.connect(dbname)
    except sqlite3.Error as e:
        print(e)
        sys.exit(-1)

    if verbose:
        sys.stderr.write("Connected to database: {}\n".format(sqlite3.version))

    return conn



def blast_to_subsys(blastf, sqlf, outf, verbose):
    """
    Convert the tab separated blast results to subsys
    :param blastf: blast/diamond output file in m8 format
    :param sqlf: sql database
    :param outf: output file
    :param verbose: more output
    :return:
    """

    """
    Table definitions:
    
    CREATE TABLE proteins(
        "patric_id" TEXT,
        "refseq_locus_tag" TEXT,
        "product" TEXT,
        "role_name" TEXT,
        "superclass" TEXT,
        "class" TEXT,
        "subclass" TEXT,
        "subsystem_name" TEXT
    );
    CREATE TABLE md5sum(
        "md5sum" TEXT,
        "patric_id" TEXT
    );

    """


    if verbose:
        message("Connecting to db", "GREEN")

    db = connect_to_db(sqlf)
    cur = db.cursor()
    if not db:
        message(f"FATAL: could not connect to database{sqlf}", "RED")
        sys.exit(-1)

    q = "select role_name,  superclass, class, subclass, subsystem_name from proteins " + \
        "inner join md5sum on proteins.patric_id = md5sum.patric_id where md5sum.md5sum = ?"

    seen = set()

    for br in stream_blast_results(blastf):
        if br.query in seen:
            continue
        seen.add(br.query)
        cur.execute(q, [br.db])
        if br.query not in proteins:
            proteins[br.query] = set()
        for r in cur.fetchall():
            """
                role = r['role_name']
                proteins[br.query].add(role)
                ss = "\t".join(r)
                if role not in roles:
                    roles[role] = set()
                roles[role].add(ss)
            """
            r.insert(0, br.query)
            print("\t".join(map(str, r)))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-f', help='blast output file', required=True)
    parser.add_argument('-d', help='database [Default: %(default)s]', default="/home/redwards/protein_functions_subsystems.sql")
    parser.add_argument('-o', help='output file to write', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    blast_to_subsys(args.f, args.d, args.o, args.v)