"""
Create an index for a fastq file with a dict of sequence ID and location in the file. This allows us
to read the fastq file and rapidly extract sequences.

We write this to a file, and provide  a method to access an index from disk
"""

import os
import sys
import argparse
import sqlite3
import gzip


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


def create_index(fqfile, indexfile, overwrite=False, verbose=False):
    """
    Index the fastq file given by fastq_file and create the index in indexfile
    :param fqfile: The fastq file to index
    :param indexfile: Where to store the output
    :param overwrite: overwrite the database file if it exists
    :param verbose: More output
    :return:
    """

    if os.path.exists(indexfile) and not overwrite:
        sys.stderr.write("Sorry, {} already exists. Please set the overwrite flag\n".format(indexfile))
        sys.exit(-1)

    conn = connect_to_db(indexfile, verbose=verbose)

    conn.execute("CREATE TABLE seqloc (id TEXT, posn INTEGER)")
    conn.commit()

    if fqfile.endswith('.gz'):
        qin = gzip.open(fqfile, 'rt')
    else:
        qin = open(fqfile, 'r')

    while qin:
        posn = qin.tell()
        header = qin.readline().split(" ")[0].replace('@', '', 1)
        conn.execute("INSERT INTO seqloc (id, posn) VALUES (header, posn)")
        seq = qin.readline()
        seq = qin.readline()
        seq = qin.readline()

    disconnect(conn, verbose=verbose)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Index a fastq file")
    parser.add_argument('-f', help='fastq file to index', required=True)
    parser.add_argument('-c', help='index file to create', required=True)
    parser.add_argument('-w', help='overwrite the file if it exists', action="store_true")
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    create_index(args.f, args.c, args.w, args.v)