import os
import sys

__author__ = 'Rob Edwards'

import sqlite3
import argparse


def get_cursor(database_file):
    """
    Open the database and get the cursor
    :param database_file:
    :type database_file:
    :return:
    :rtype:
    """

    sql = sqlite3.connect(database_file)
    return sql.cursor()

def count_study_types(database_file):
    """
    Count the study types in the database

    :param database_file:
    :type database_file:
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



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='get the study types')
    parser.add_argument('-d', help='SQLLite databasae')
    args = parser.parse_args()

    count_study_types(args.d)