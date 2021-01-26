"""
Convert the NCBI taxonomy to an SQLite database for faster access
"""

import os
import sys
import argparse
import sqlite3
import json

from config import get_db_dir
import gzip

defaultdir = get_db_dir()

def connect_to_db(dbname, verbose=False):
    """
    Connect to the database
    :param dbname: the database file name
    :param verbose: print addtional output
    :return: the database connection
    """

    try:
        if verbose:
            sys.stderr.write("Connecting to {}\n".format(os.path.join(defaultdir, dbname)))
        conn = sqlite3.connect(os.path.join(defaultdir, dbname))
    except sqlite3.Error as e:
        sys.stderr.write("ERROR Creating database: {}\n".format(os.path.join(defaultdir, dbname)))
        sys.stderr.write(e)
        sys.exit(-1)

    if verbose:
        sys.stderr.write("Connected to database: {}\n".format(sqlite3.version))

    return conn


def disconnect(conn, verbose=False):
    """
    Disconnect the database and ensure we've saved all changes
    :param conn: the database connection
    :param verbose: print addtional output
    :return:
    """

    if conn:
        conn.commit()
        conn.close()
    elif verbose:
        sys.stderr.write("There was no database connection!\n")

def create_load(conn, datadir, verbose=False):
    """
    Create the databases and load the data.

    Note: Initially I wrote a generic version of this, but you need to have the correct
    number of ?'s for this. We probably could generalize some of this, especially those
    with two columns, but this works well.

    :param conn: the database connection
    :param datadir: the databae directory
    """


    ## The NODES table
    dbfile = os.path.join(datadir, "nodes.dmp")
    if verbose:
        sys.stderr.write("loading NODES table: {}\n".format(dbfile))
    if not os.path.exists(dbfile):
        sys.stderr.write("ERROR: {} does not exist\n".format(dbfile))
        sys.exit(-1)
    conn.execute("CREATE TABLE nodes (tax_id INTEGER PRIMARY KEY, parent INTEGER, rank TEXT, embl_code TEXT, division_id INTEGER, inherited_div INTEGER, genetic_code INTEGER, inherited_genetic_code INTEGER, mitochondrial_genetic_code INTEGER, inherited_mito_gc INTEGER, genbank_hidden INTEGER, hidden_subtree INTEGER, comments)")
    conn.commit()
    with open(dbfile, "r") as f:
        for l in f:
            p = l.strip().rstrip("|").split('\t|')
            p = [x.strip() for x in p]
            try:
                conn.execute("INSERT INTO nodes VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", p)
            except sqlite3.OperationalError as e:
                sys.stderr.write("{}".format(e))
                sys.stderr.write("\nWhile insert on: {}\n".format(p))
                sys.exit()
    conn.commit()


    # the NAMES table
    dbfile = os.path.join(datadir, "names.dmp")
    if verbose:
        sys.stderr.write("loading NAMES table: {}\n".format(dbfile))
    if not os.path.exists(dbfile):
        sys.stderr.write("ERROR: {} does not exist\n".format(dbfile))
        sys.exit(-1)
    conn.execute("CREATE TABLE names (tax_id INTEGER, name TEXT, unique_name TEXT, name_class TEXT)")
    conn.commit()
    with open(dbfile, "r") as f:
        for l in f:
            p = l.strip().rstrip("|").split('\t|')
            p = [x.strip() for x in p]
            try:
                conn.execute("INSERT INTO names VALUES (?, ?, ?, ?)", p)
            except sqlite3.OperationalError as e:
                sys.stderr.write("{}".format(e))
                sys.stderr.write("\nWhile insert on: {}\n".format(p))
                sys.exit()
    conn.commit()

    # the DIVISION table
    dbfile = os.path.join(datadir, "division.dmp")
    if verbose:
        sys.stderr.write("loading DIVISION table: {}\n".format(dbfile))
    if not os.path.exists(dbfile):
        sys.stderr.write("ERROR: {} does not exist\n".format(dbfile))
        sys.exit(-1)
    conn.execute("CREATE TABLE division (division_id INTEGER PRIMARY KEY, division_code TEXT, division_name TEXT, comments)")
    conn.commit()
    with open(dbfile, "r") as f:
        for l in f:
            p = l.strip().rstrip("|").split('\t|')
            p = [x.strip() for x in p]
            try:
                conn.execute("INSERT INTO division VALUES (?, ?, ?, ?)", p)
            except sqlite3.OperationalError as e:
                sys.stderr.write("{}".format(e))
                sys.stderr.write("\nWhile insert on: {}\n".format(p))
                sys.exit()
    conn.commit()

    # the GENETIC CODE table
    dbfile = os.path.join(datadir, "gencode.dmp")
    if verbose:
        sys.stderr.write("loading GENETIC CODE table: {}\n".format(dbfile))
    if not os.path.exists(dbfile):
        sys.stderr.write("ERROR: {} does not exist\n".format(dbfile))
        sys.exit(-1)
    conn.execute("CREATE TABLE gencode (genetic_code INTEGER, abbreviation , name TEXT, cde TEXT, starts TEXT)")
    conn.commit()
    with open(dbfile, "r") as f:
        for l in f:
            p = l.strip().rstrip("|").split('\t|')
            p = [x.strip() for x in p]
            try:
                conn.execute("INSERT INTO gencode VALUES (?, ?, ?, ?, ?)", p)
            except sqlite3.OperationalError as e:
                sys.stderr.write("{}".format(e))
                sys.stderr.write("\nWhile insert on: {}\n".format(p))
                sys.exit()
    conn.commit()

    # the MERGED database
    dbfile = os.path.join(datadir, "merged.dmp")
    if verbose:
        sys.stderr.write("loading MERGED table: {}\n".format(dbfile))
    if not os.path.exists(dbfile):
        sys.stderr.write("ERROR: {} does not exist\n".format(dbfile))
        sys.exit(-1)
    conn.execute("CREATE TABLE merged (old_tax_id INTEGER, new_tax_id INTEGER)")
    conn.commit()
    with open(dbfile, "r") as f:
        for l in f:
            p = l.strip().rstrip("|").split('\t|')
            p = [x.strip() for x in p]
            try:
                conn.execute("INSERT INTO merged VALUES (?, ?)", p)
            except sqlite3.OperationalError as e:
                sys.stderr.write("{}".format(e))
                sys.stderr.write("\nWhile insert on: {}\n".format(p))
                sys.exit()
    conn.commit()

    # The deleted nodes
    dbfile = os.path.join(datadir, "delnodes.dmp")
    if verbose:
        sys.stderr.write("loading DELETED table: {}\n".format(dbfile))
    if not os.path.exists(dbfile):
        sys.stderr.write("ERROR: {} does not exist\n".format(dbfile))
        sys.exit(-1)
    conn.execute("CREATE TABLE deleted (tax_id INTEGER)")
    conn.commit()
    with open(dbfile, "r") as f:
        for l in f:
            l = l.replace("|", "").strip()
            try:
                conn.execute("INSERT INTO deleted VALUES (?)", [l])
            except sqlite3.OperationalError as e:
                sys.stderr.write("{}".format(e))
                sys.stderr.write("\nWhile insert on: {}\n".format(l))
                sys.exit()
    conn.commit()

    return conn

def not_used(conn, datadir, verbose=False):
    """
    Unreachable code!
    
    The following two sections are deprecated because NCBI no longer supports gi_taxid tables directly
    
    However, there is the new accession2taxid folder with similar data and we should update this code
    with that folder.

    This is here as legacy code
    """

    # the gi_taxid
    dbfile = os.path.join(datadir, "gi_taxid_nucl.dmp.gz")
    if verbose:
        sys.stderr.write("loading GI_TAXID (NUCL) table: {}\n".format(dbfile))
    if not os.path.exists(dbfile):
        sys.stderr.write("ERROR: {} does not exist\n".format(dbfile))
        sys.exit(-1)
    conn.execute("CREATE TABLE gi_taxid_nucl (gi INTEGER PRIMARY KEY, tax_id INTEGER)")
    conn.commit()
    with gzip.open(dbfile, "rt") as f:
        for l in f:
            p = l.strip().split('\t')
            p = [x.strip() for x in p]
            try:
                conn.execute("INSERT INTO gi_taxid_nucl VALUES (?, ?)", p)
            except sqlite3.OperationalError as e:
                sys.stderr.write("{}".format(e))
                sys.stderr.write("\nWhile insert on: {}\n".format(p))
                sys.exit()
    conn.commit()

    # the gi_taxid
    dbfile = os.path.join(datadir, "gi_taxid_prot.dmp.gz")
    if verbose:
        sys.stderr.write("loading GI_TAXID (PROT) table: {}\n".format(dbfile))
    if not os.path.exists(dbfile):
        sys.stderr.write("ERROR: {} does not exist\n".format(dbfile))
        sys.exit(-1)
    conn.execute("CREATE TABLE gi_taxid_prot (gi INTEGER PRIMARY KEY, tax_id INTEGER)")
    conn.commit()
    with gzip.open(dbfile, "rt") as f:
        for l in f:
            p = l.strip().split('\t')
            p = [x.strip() for x in p]
            try:
                conn.execute("INSERT INTO gi_taxid_prot VALUES (?, ?)", p)
            except sqlite3.OperationalError as e:
                sys.stderr.write("{}".format(e))
                sys.stderr.write("\nWhile insert on: {}\n".format(p))
                sys.exit()
    conn.commit()


    return conn

def create_indices(conn, verbose=False):
    """
    Create some useful indices. Note that the PRIMARY KEY columns are indexed by default!
    :param conn: The database connection
    :param verbose: print addtional output
    :return:
    """
    
    if verbose:
        sys.stderr.write("Creating indices\n")

    tables = {
        "nodes": {"tidparentrank" : ["tax_id", "parent", "rank"]},
        "names": {
            "tidname" : ["tax_id", "name"],
            "tiduniname" : ["tax_id", "unique_name"],
            "tidnameuniname" : ["tax_id", "name", "unique_name"]
        },
        "division": {"divname" : ["division_id", "division_name"]},
        "merged" : {"oldnewidx" : ["old_tax_id", "new_tax_id"]}
    }

    for t in tables:
        for idx in tables[t]:
            conn.execute("CREATE INDEX {ix} ON {tn} ({cn})".format(ix=idx, tn=t, cn=", ".join(tables[t][idx])))
    conn.commit()

    return conn


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert NCBI taxonomy to SQLite3 database")
    parser.add_argument('-d', help='directory with tax data', required=True)
    parser.add_argument('-s', help='sqlite file to write (default=taxonomy.sqlite3)',default="taxonomy.sqlite3")
    parser.add_argument('-v', help='verbose output', action="store_true")
    parser.add_argument('-o', help='overwrite any existing databases (otherwise error out)', action="store_true")
    args = parser.parse_args()

    if os.path.exists(args.s):
        if args.o:
            os.remove(args.s)
        else:
            sys.stderr.write("{} already exists so we won't overwrite it. Use -o to force overwrite\n".format(args.s))
            sys.exit(-1)

    conn = connect_to_db(args.s, args.v)
    conn = create_load(conn, args.d, args.v)
    conn = create_indices(conn, args.v)
    disconnect(conn, args.v)
