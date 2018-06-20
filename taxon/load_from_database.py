"""
Lazy loading of taxonomy data from the sqlite3 database
"""

import os
import sys
import sqlite3
import argparse

from .taxonomy import TaxonNode, TaxonName, TaxonDivision

data = {"node": {}, "name": {}, "division": {}}
default_database = "/raid60/usr/data/NCBI/taxonomy/current/taxonomy.sqlite3"

def connect_to_db(dbname, verbose=False):
    """
    Connect to the database
    :param dbname: the database file name
    :param verbose: print addtional output
    :return: the database connection
    """

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

def get_taxonomy(taxid, conn):
    """
    Retrieve a TaxonNode object for a given taxonomy ID
    :param taxid: the taxonomy id
    :param conn: the database connection
    :return: a TaxonNode and TaxonName object
    """

    global data
    cur = conn.cursor()

    if taxid in data['node']:
        return data['node'][taxid]
    else:
        sys.stderr.write("Trying to find {}\n".format(taxid))
        cur.execute("select * from nodes where tax_id = ?", [taxid])
        p = cur.fetchone()
        t = TaxonNode(*p)
        data['node'][taxid] = t


        cur.execute("select * from names where tax_id = ?", [taxid])
        n = TaxonName(taxid)
        for p in cur.fetchall():
            if p[2]:
                n.unique = p[2]
            else:
                n.set_name(p[3], p[1])
        
        return t, n



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get data from the database")
    parser.add_argument('-s', help='sqlite file to write',required=True)
    parser.add_argument('-t', help='taxid to get')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    conn = connect_to_db(args.s)
    t,n = get_taxonomy(args.t, conn)
    print("{}".format(t.rank))
    print("{}: scientific: {} common: {} blast: {}\n".format(
        n.taxid, n.scientific_name, n.common_name, n.blast_name))
    disconnect(conn)
