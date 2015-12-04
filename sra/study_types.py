import os
import sys

__author__ = 'Rob Edwards'

import sqlite3
import argparse


def get_cursor(database_file):
    """
    Open the database and get the cursor
    :param database_file:The sql lite database
    :type database_file:str
    :return:
    :rtype:
    """

    sql = sqlite3.connect(database_file)
    return sql.cursor()

def count_study_types(database_file):
    """
    Count the study types in the database

    :param database_file:The sql lite database
    :type database_file:str
    :return:
    :rtype:
    """

    sql = get_cursor(database_file)
    count = {}
    for row in sql.execute('select study_type from study'):
        count[row[0]] = count.get(row[0], 0)+1

    print("Study type\tInstances")
    for t in count:
        print("{}\t{}".format(t, count[t]))


def count_metagenome_types(database_file):
    """
    Count the library strategies with the different metagenome datasets
    :param database_file: The sql lite database
    :type database_file: str
    :return:
    :rtype:
    """

    sql = get_cursor(database_file)
    count = {}
    print("Amplicon Studies")
    for row in sql.execute(
            'select study_type, count(1) from study where study_accession in (select study_accession from experiment where library_strategy = "AMPLICON") group by study_type;'):
        print("\t".join(row))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='get the study types')
    parser.add_argument('-d', help='SQLLite databasae')
    args = parser.parse_args()

    count_study_types(args.d)
    count_metagenome_types(args.d)