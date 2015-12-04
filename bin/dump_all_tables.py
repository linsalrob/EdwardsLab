import argparse
import sqlite3
import pandas as pd


def to_csv(filename):
    db = sqlite3.connect(filename)
    cursor = db.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cursor.fetchall()
    for table_name in tables:
        table_name = table_name[0]
        table = pd.read_sql_query("SELECT * from %s" % table_name, db)
        table.to_csv(table_name + '.csv', index_label='index')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Dump the contents of an SQL file to CSV. This was taken from http://stackoverflow.com/questions/305378/get-list-of-tables-db-schema-dump-etc-in-sqlite-databases')
    parser.add_argument('-d', help='SQLlite database file', required=True)
    args = parser.parse_args()
    to_csv(args)
