"""
Parse the model seed biochemistry file and load it into an SQL lite table
"""

import os
import sys
import argparse
import sqlite3
import json

from roblib import bcolors

__author__ = 'Rob Edwards'

ModelSeedLocation = os.environ['ModelSEEDDatabase']

if not ModelSeedLocation:
    sys.stderr.write(f"{bcolors.RED}FATAL{bcolors.ENDC} The {bcolors.BLUE}ModelSEEDDatabase{bcolors.ENDC} environment variable is not set\n")
    sys.exit(255)

if not os.path.exists(ModelSeedLocation):
    sys.stderr.write(f"{bcolors.RED}FATAL{bcolors.ENDC} The {bcolors.BLUE}ModelSEEDDatabase{bcolors.ENDC} environment variable is set but the path does not exist\n")
    sys.exit(255)


def createdb(db, verbose=False):
    """
    Create the database connection
    :param db:
    :return:
    """

    try:
        conn = sqlite3.connect(db)
    except sqlite3.Error as e:
        sys.stderr.write("ERROR Creating database: {}\n".format(db))
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


def create_tables(conn, verbose=False):
    """
    Create a table for a given genomeid
    :param conn: the database connection
    :param verbose: more output
    :return: the database connection
    """

    # compounds
    rows = 'charge INTEGER, name TEXT, id TEXT, smiles TEXT, abstract_compound TEXT, is_obsolete INTEGER, '
    rows += 'linked_compound TEXT, pkb TEXT, source TEXT, abbreviation TEXT, comprised_of TEXT, pka TEXT, '
    rows += 'deltagerr REAL, formula TEXT, is_core TEXT, is_cofactor TEXT, deltag REAL, mass INTEGER, inchikey TEXT, aliases TEXT'


    cl = "CREATE TABLE compounds ({rows})".format(rows=rows)
    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}CREATING compounds table{bcolors.ENDC}\n")
        sys.stderr.write("\t{}\n".format(cl))
    conn.execute(cl)






    conn.commit()

    return conn


def create_database(dbfile, verbose=False):
    """
    Create the SQL Lite database at the location of dbfile
    :param db: the SQL lite database to create
    :param verbose: More output
    :return: nothing
    """

    c = createdb(dbfile, verbose)

    # load the compounds
    cpds = json.load(open(os.path.join(ModelSeedLocation, "compounds.json"), 'r'))
    for c in cpds:







if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='input file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()