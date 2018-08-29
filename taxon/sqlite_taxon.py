"""
Convert the NCBI taxonomy to an SQLite database for faster access
"""

import os
import sys
import argparse
import sqlite3
import json

from config import get_db_dir

defaultdir = get_db_dir()

def connect_to_db(dbname, verbose=False):
    """
    Connect to the database
    :param dbname: the database file name
    :param verbose: print addtional output
    :return: the database connection
    """

    try:
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

def define_tables(conn, verbose=False):
    """
    Define the database
    :param conn: the database connection
    :param verbose: print addtional output
    :return: the database connection
    """
    """
    Tables in the database (note that the tables have the same name as the files on purpose!):
    
    Note: You can find more information about these files at ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_readme.txt
    
    nodes:
            tax_id                                  -- node id in GenBank taxonomy database
            parent tax_id                           -- parent node id in GenBank taxonomy database
            rank                                    -- rank of this node (superkingdom, kingdom, ...)
            embl code                               -- locus-name prefix; not unique
            division id                             -- see division.dmp file
            inherited div flag  (1 or 0)            -- 1 if node inherits division from parent
            genetic code id                         -- see gencode.dmp file
            inherited GC  flag  (1 or 0)            -- 1 if node inherits genetic code from parent
            mitochondrial genetic code id           -- see gencode.dmp file
            inherited MGC flag  (1 or 0)            -- 1 if node inherits mitochondrial gencode from parent
            GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
            hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
            comments                                -- free-text comments and citations
        
    
    names:
            tax_id                                  -- the id of node associated with this name
            name_txt                                -- name itself
            unique name                             -- the unique variant of this name if name not unique
            name class                              -- (synonym, common name, ...)
    
    division:
            division id                             -- taxonomy database division id
            division cde                            -- GenBank division code (three characters)
            division name                           -- e.g. BCT, PLN, VRT, MAM, PRI...
            comments
    
    gencode:
            genetic code id                         -- GenBank genetic code id
	        abbreviation                            -- genetic code name abbreviation
	        name                                    -- genetic code name
	        cde                                     -- translation table for this genetic code
	        starts                                  -- start codons for this genetic code
        
    """


    tables = {
        "nodes" : {
            "tax_id"                     : "INTEGER",
            "parent"                     : "INTEGER",
            "rank"                       : "TEXT",
            "embl_code"                  : "TEXT",
            "division_id"                : "INTEGER",
            "inherited_div"              : "INTEGER",
            "genetic_code"               : "INTEGER",
            "inherited_genetic_code"     : "INTEGER",
            "mitochondrial_genetic_code" : "INTEGER",
            "inherited_mito_gc"          : "INTEGER",
            "genbank_hidden"             : "INTEGER",
            "hidden_subtree"             : "INTEGER",
            "comments" : ""
        },
        "names" : {
            "tax_id"                     : "INTEGER",
            "name"                       : "TEXT",
            "unique_name"                : "TEXT",
            "name_class"                 : "TEXT"
        },
        "division" : {
            "division_id"                : "INTEGER",
            "division_code"              : "TEXT",
            "division_name"              : "TEXT",
            "comments"                   : ""
        },
        "gencode" : {
            "genetic_code"               : "INTEGER",
            "abbreviation"               : "",
            "name"                       : "TEXT",
            "cde"                        : "TEXT",
            "starts"                     : "TEXT"
        }
    }

    # note that tax_id is not unique in names
    primary_keys = {"nodes" : "tax_id", "division" : "division_id", "gencode" : "genetic_code", "names" : ""}
    for t in tables:
        createl = []
        for c in tables[t]:
            if primary_keys[t] == c:
                createl.append("{} {} PRIMARY KEY".format(c, tables[t][c]))
            else:
                createl.append("{} {}".format(c, tables[t][c]))
        if verbose:
            sys.stderr.write("Going to create table with command: ")
            sys.stderr.write("CREATE TABLE {tn} ({rows})".format(tn=t, rows=", ".join(createl)))
            sys.stderr.write("\n")
        conn.execute("CREATE TABLE {tn} ({rows})".format(tn=t, rows=", ".join(createl)))
    conn.commit()

    return conn


def load_tables(conn, datadir, verbose=False):
    """
    Parse the files and load the tables
    :param conn: the database connection
    :param datadir:  the data directory with the files
    :param verbose: print addtional output
    :return: the database connection
    """

    tables = {
        "nodesXX": ["tax_id", "parent", "rank", "embl_code", "division_id", "inherited_div", "genetic_code",
                      "inherited_genetic_code", "mitochondrial_genetic_code", "inherited_mito_gc", "genbank_hidden",
                      "hidden_subtree", "comments"],
        "names": ["tax_id", "name", "unique_name", "name_class"],
        "division": ["division_id", "division_code", "division_name", "comments"],
        "gencode": ["genetic_code", "abbreviation", "name", "cde", "starts"]
    }


    for t in tables:
        if not os.path.exists(os.path.join(datadir, "{}.dmp".format(t))):
            sys.stderr.write("ERROR: {}.dmp not found in {}\n".format(t, datadir))
            continue
        vcount = 0
        with open(os.path.join(datadir, "{}.dmp".format(t)), 'r') as f:
            for l in f:
                pori = l.strip().rstrip("|").split('\t|')
                p=[]
                for i,j in enumerate(pori):
                    x=j.replace('"', '')
                    x=x.strip()
                    p.append(x)
                sys.stderr.write("{}\n".format(p))
                # marvel in this json.dumps!
                maxvcount  = 100
                if verbose:
                    vcount += 1
                    if vcount < maxvcount:
                        sys.stderr.write("INSERT INTO {tn} ({cn}) VALUES ({val})\n".format(
                        tn=t, cn=", ".join(tables[t]), val=json.dumps(p).replace('[', '', 1).replace(']', '', 1)))
                    if vcount == maxvcount:
                        sys.stderr.write("Only first {mv} inputs printed\n\n".format(mv=maxvcount))
                try:
                    conn.execute("INSERT INTO {tn} ({cn}) VALUES ({val})".format(
                        tn=t, cn=", ".join(tables[t]), val=json.dumps(p).replace('[', '', 1).replace(']', '', 1)
                    ))
                except sqlite3.OperationalError as e:
                    sys.stderr.write("{}".format(e))
                    sys.stderr.write("\nWhile performing:\n")
                    sys.stderr.write("INSERT INTO {tn} ({cn}) VALUES ({val})\n".format(
                        tn=t, cn=", ".join(tables[t]), val=json.dumps(p).replace('[', '', 1).replace(']', '', 1)))
                    sys.exit()

        conn.commit()
    return conn


def create_load(conn, datadir, verbose=False):
    """
    Create the databases and load the data. This is not generic like the
    other methods, but it allows for quotes etc in the names. grr.
    :param conn: the database connection
    :param datadir: the databae directory
    """


    ## The NODES table
    if verbose:
        sys.stderr.write("loading nodes table\n")
    conn.execute("CREATE TABLE nodes (tax_id INTEGER PRIMARY KEY, parent INTEGER, rank TEXT, embl_code TEXT, division_id INTEGER, inherited_div INTEGER, genetic_code INTEGER, inherited_genetic_code INTEGER, mitochondrial_genetic_code INTEGER, inherited_mito_gc INTEGER, genbank_hidden INTEGER, hidden_subtree INTEGER, comments)")
    conn.commit()
    with open(os.path.join(datadir, "nodes.dmp"), "r") as f:
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
    if verbose:
        sys.stderr.write("loading names table\n")
    conn.execute("CREATE TABLE names (tax_id INTEGER, name TEXT, unique_name TEXT, name_class TEXT)")
    conn.commit()
    with open(os.path.join(datadir, "names.dmp"), "r") as f:
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
    if verbose:
        sys.stderr.write("loading division table\n")
    conn.execute("CREATE TABLE division (division_id INTEGER PRIMARY KEY, division_code TEXT, division_name TEXT, comments)")
    conn.commit()
    with open(os.path.join(datadir, "division.dmp"), "r") as f:
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
    if verbose:
        sys.stderr.write("loading gencode table\n")
    conn.execute("CREATE TABLE gencode (genetic_code INTEGER, abbreviation , name TEXT, cde TEXT, starts TEXT)")
    conn.commit()
    with open(os.path.join(datadir, "gencode.dmp"), "r") as f:
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
    if verbose:
        sys.stderr.write("loading merged table\n")
    conn.execute("CREATE TABLE merged (old_tax_id INTEGER, new_tax_id INTEGER)")
    conn.commit()
    with open(os.path.join(datadir, "merged.dmp"), "r") as f:
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
    parser.add_argument('-s', help='sqlite file to write',required=True)
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
    # conn = define_tables(conn, args.v)
    # conn = load_tables(conn, args.d, args.v)
    conn = create_load(conn, args.d, args.v)
    conn = create_indices(conn, args.v)
    disconnect(conn, args.v)
