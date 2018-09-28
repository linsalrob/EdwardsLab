"""

"""

import os
import sys
import argparse
import sqlite3
import gzip

__author__ = 'Rob Edwards'

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


def create_table(conn, verbose=False):
    """
    Create a table for a given genomeid
    :param conn: the database connection
    :param verbose: more output
    :return: the database connection
    """

    rows =  'query TEXT, subject TEXT, identity REAL, len INTEGER, mismatch INTEGER,'
    rows += ' gap_open INTEGER, query_start INTEGER, query_end INTEGER, subject_start INTEGER,'
    rows += ' subject_end INTEGER, e_value REAL, bitscore REAL'

    cl = "CREATE TABLE rapsearch ({rows})".format(rows=rows)
    if verbose:
        sys.stderr.write("{}\n".format(cl))

    conn.execute(cl)
    conn.execute("CREATE INDEX rapsearchQIndex on rapsearch(query)")
    conn.execute("CREATE INDEX rapsearchSIndex on rapsearch(subject)")

    conn.commit()

    return conn

def load_data(filename, conn, verbose=False):
    """
    Load the data fromt the m8 format file
    :param filename: the file for the m8 format data
    :param conn: the database connection
    :param verbose: more output
    :return:
    """

    if verbose:
        sys.stderr.write("Parsing {}\n".format(filename))

    if filename.endswith('.gz'):
        qin = gzip.open(filename, 'rt')
    else:
        qin = open(filename, 'r')

    for l in qin:
        if l.startswith('#'):
            continue
        p = l.strip().split("\t")
        """
        for i  in [2, 10, 11]:
            p[i] = float(p[i])
        for i in [3, 4, 5, 6, 7, 8, 9]:
            p[i] = int(p[i])
        """

        c = "INSERT INTO rapsearch (query, subject, identity, len, mismatch, gap_open, query_start, query_end, subject_start, subject_end, e_value, bitscore)"
        c += "values ('{val}')".format(val="', '".join(p))

        try:
            conn.execute(c)
        except sqlite3.OperationalError as e:
            sys.stderr.write("{}".format(e))
            sys.stderr.write("\nWhile performing:\n{}\n".format(c))
            sys.exit(-1)

        conn.commit()

    qin.close()

    return conn

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-d', help='directory of rapsearch files', required=True)
    parser.add_argument('-s', help='sqlite database to write to')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    conn = createdb(args.s, args.v)
    create_table(conn, args.v)

    for f in os.listdir(args.d):
        conn = load_data(os.path.join(args.d, f), conn, args.v)
    disconnect(conn)