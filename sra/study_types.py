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
        print("\t".join(map(str, row)))



def get_metagenome_run_ids(database_file, amplicons=False):
    """
    Get all the run id's for the metagenome sequences, If amplicons is false we get the non-amplicon sequences. If
    amplicons is true we only get the amplicon sequences.

    :param database_file: The sql lite file
    :type database_file: str
    :param amplicons: Whether to get amplicon sequences or not
    :type amplicons: bool
    :return:
    :rtype:
    """

    sql = get_cursor(database_file)
    # select experiment_accession from study left join experiment on
    # study.study_accession = experiment.study_accession where experiment.library_strategy = 'FL-cDNA' and study_type NOT LIKE '%their%';

    # first build the query to get the experiment accession ids
    query = 'select experiment_accession from study left join experiment on study.study_accession = '
    query += 'experiment.study_accession where experiment.library_strategy'

    if amplicons:
        query += ' = "AMPLICON"'
    else:
        query += ' != "AMPLICON"'
    query += ' and study_type = "Metagenomics"'

    # now wrap that to get the accession ids using a subselect
    query = 'select run_accession from run where experiment_accession in (' + query + ');'

    sys.stderr.write("Executing query\n" + query + "\n")

    for row in sql.execute(query):
        print(str(row[0]))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='get the study types')
    parser.add_argument('-d', help='SQLLite databasae')
    args = parser.parse_args()

    # count_study_types(args.d)
    # count_metagenome_types(args.d)
    print("NON AMPLICON RUNS")
    get_metagenome_run_ids(args.d, False)
    print("\nAMPLICON RUNS")
    get_metagenome_run_ids(args.d, True)
