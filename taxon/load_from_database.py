"""
Lazy loading of taxonomy data from the sqlite3 database
"""

import os
import sys
import sqlite3
import argparse

from .taxonomy import TaxonNode, TaxonName, TaxonDivision
from .Error import EntryNotInDatabaseError
from roblib import bcolors
import time

data = {"node": {}, "name": {}, "division": {}}
default_database = "/raid60/usr/data/NCBI/taxonomy/current/taxonomy.sqlite3"
conn = None
#default_database = "/data/ncbi/taxonomy/20180620/taxonomy.sqlite"

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



def get_taxonomy_db():
    """
    Connect to the default SQLite3 taxonomy database
    """
    
    if os.path.exists(default_database):
        return connect_to_db(default_database)
    else:
        sys.stderr.write("The default database ({}) does not exist. Please create a connection\n".format(default_database))
        sys.exit(-1)

def get_taxid_for_name(name, conn, verbose=False):
    """
    Retrieve the taxnomy ID for a scientific name
    :param name: The name to look up
    :param conn: The database connection
    :param verbose: more output
    :return: The taxonomy ID
    """

    cur = conn.cursor()
    cur.execute("select tax_id from names where name = ?", [name])

    newid = cur.fetchone()
    if newid and newid[0]:
        return newid[0]
    else:
        sys.stderr.write(f"{bcolors.FAIL}ERROR:{bcolors.ENDC}: {name} is not in the database\n")
        return None


def get_taxonomy(taxid, conn, verbose=False):
    """
    Retrieve a TaxonNode object for a given taxonomy ID
    :param taxid: the taxonomy id
    :param conn: the database connection
    :param verbose: more output
    :return: a TaxonNode and TaxonName object
    """

    global data
    cur = conn.cursor()
    if taxid in data['node']:
        return data['node'][taxid], data['name'][taxid]

    cur.execute("select * from nodes where tax_id = ?", [taxid])
    p = cur.fetchone()
    if not p:
        # check the merged database
        cur.execute("select new_tax_id from merged where old_tax_id = ?", [taxid])
        newid = cur.fetchone()
        if newid and newid[0]:
            cur.execute("select * from nodes where tax_id = ?", [newid[0]])
            p = cur.fetchone()
        else:
            raise EntryNotInDatabaseError(f"ERROR: {taxid} is not in the database and not merged\n")

    t = TaxonNode(*p)
    data['node'][taxid] = t


    cur.execute("select * from names where tax_id = ?", [taxid])
    n = TaxonName(taxid)
    for p in cur.fetchall():
        if p[2]:
            n.unique = p[2]
        n.set_name(p[3], p[1])
    data['name'][taxid] = n
    return t, n


def gi_to_taxonomy(gi, conn, protein=False, verbose=False):
    """
    Convert an NCBI gi to a taxonomy object
    :param gi: The NCBI gi number
    :param conn: the database connection
    :param protein: Whether the object refers to protein (True) or DNA (False). Default=DNA
    :param verbose: More output
    :return: A taxonomy object
    """

    global data
    cur = conn.cursor()
    if gi in data['gi2tax']:
        taxid = data['gi2tax']
        return data['node'][taxid], data['name'][taxid]

    db = "gi_taxid_nucl"
    if protein:
        db = "gi_taxid_prot"
    sqlexe="select tax_id from {} where gi = ?".format(db)
    cur.execute(sqlexe, [gi])
    p = cur.fetchone()[0]
    data['gi2tax'][gi] = p
    if verbose:
        sys.stderr.write("GI: {} Taxonomy: {}\n".format(gi, p))
    return get_taxonomy(p, conn)

def all_ids_complete(conn, protein=False, verbose=False):
    """
    Get all the available IDs in the database
    :param conn: the database connection
    :param protein: Whether the object refers to protein (True) or DNA (False). Default=DNA
    :param verbose: More output
    :param verbose: More output
    :return: A list of taxonomy ids
    """
    global data
    cur = conn.cursor()

    cur.execute("select * from nodes")
    sys.stderr.write(f"{bcolors.YELLOW}Collecting all the data. Please stand by.\n{bcolors.ENDC}")
    sys.stderr.write(f"{bcolors.RED}Warning, this will take a long time!!.\n{bcolors.ENDC}")
    for p in cur.fetchall():
        t = TaxonNode(*p)
        data['node'][p[0]] = t
        cur.execute("select * from names where tax_id = ?", [p[0]])
        n = TaxonName(p[0])
        for r in cur.fetchall():
            if r[2]:
                n.unique = r[2]
            n.set_name(r[3], r[1])
        data['name'][p[0]] = n
    sys.stderr.write(f"{bcolors.GREEN}Done.\n{bcolors.ENDC}")
    return t, n

def all_ids(conn, verbose=False):
    """
    Just return a list of all taxonomy IDs
    :param conn: database connection
    :param verbose: more output
    :return: a list of just the taxonomy IDs
    """

    cur = conn.cursor()
    sys.stderr.write(f"{bcolors.YELLOW}Collecting all the taxon data. Please stand by.\n{bcolors.ENDC}")
    sys.stderr.write(f"{bcolors.RED}Warning, this will take a long time!!.\n{bcolors.ENDC}")
    s = time.time()
    exc = cur.execute("select tax_id from nodes")
    sys.stderr.write(f"{bcolors.GREEN}Done in {time.time() - s} seconds!.\n{bcolors.ENDC}")
    return exc.fetchall()

def all_species_ids(conn, verbose=False):
    """
    Get all the species taxonomy IDs.

    We filter for rank == species, but only return the taxonomy ID
    :param conn: the database connection
    :param verbose: more output
    :return: all of the taxonomy IDs
    """
    cur = conn.cursor()
    exc = cur.execute("select tax_id from nodes where rank='species'")
    return exc.fetchall()


def taxonomy_hierarchy(tid, verbose=False):
    """
    Get the taxonomical hierarchy for a tax id. Yields so you can call this in a while loop
    Note we just yield the id
    :param tid: taxonomy ID
    :param verbose: More output
    """

    global data
    global conn


    while tid != 1:
        if tid == 1:
            if verbose:
                sys.stderr.write(f"{bcolors.RED}Taxonomy ID 1 found{bcolors.ENDC}\n")
            return
        if not tid:
            if verbose:
                sys.stderr.write(f"{bcolors.RED}No tid{bcolors.ENDC}\n")
            return

        if tid not in data['node']:
            try:
                get_taxonomy(tid, conn, verbose)
            except EntryNotInDatabaseError:
                if verbose:
                    sys.stderr.write(f"{bcolors.RED}{tid} is not in database. Can not continue{bcolors.ENDC}\n")
                return

        if verbose:
            sys.stderr.write(f"{bcolors.GREEN}tid: {tid} parent: {data['node'][tid].parent}{bcolors.ENDC}\n")
        yield data['node'][tid].parent
        tid = data['node'][tid].parent

def taxonomy_hierarchy_as_list(conn, tid, verbose=False):
    """
    Return the taxonomy hierarchy as a list
    :param conn: the database connect
    :param tid: the taxonomy id
    :param verbose: more output
    :return:
    """
    wanted_levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies']
    taxlist = ["", "", "", "", "", "", "", ""]

    t, n = get_taxonomy(tid, conn)
    if not t:
        if verbose:
            sys.stderr.write("No taxonomy for {}\n".format(tid))
        return taxlist

    while t.parent != 1 and t.taxid != 1:
        if t.rank in wanted_levels:
            taxlist[wanted_levels.index(t.rank)] = n.scientific_name
        t, n = get_taxonomy(t.parent, conn)
    return taxlist

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get data from the database")
    parser.add_argument('-s', help='sqlite database file',required=True)
    parser.add_argument('-t', help='taxid to get')
    parser.add_argument('-g', help='gi to get', type=int)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    conn = connect_to_db(args.s)

    if args.t:
        t,n = get_taxonomy(args.t, conn, args.v)
    elif args.g:
        t,n = gi_to_taxonomy(args.g, conn, verbose=args.v)
    else:
        sys.stderr.write("Please provide one of either -t or -g")
        sys.exit(-1)

    print("{}".format(t.rank))
    print("{}: scientific: {} common: {} blast: {}\n".format(
        n.taxid, n.scientific_name, n.common_name, n.blast_name))
    disconnect(conn)
