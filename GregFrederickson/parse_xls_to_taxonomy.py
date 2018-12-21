"""
Parse the meta_abundance.xls file from Kate (see email, Dec 6th) and merge it with the NCBI taxonomy.
"""

import os
import sys
import argparse
from openpyxl import load_workbook
from taxon import connect_to_db, get_taxonomy, get_taxid_for_name


__author__ = 'Rob Edwards'

taxonomy_str = {}

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def parse_sheet(xlsfile, verbose=True, nocolor=False):
    """
    Read the xls file and extract the data we want
    :param xlsfile: the file to read
    :param verbose: more output
    :param nocolor: no color ouput
    :return:
    """

    correct_names = {"Clostridium beijerincki": "Clostridium beijerinckii",
                     "Caldicellulosiruptor accharolyticus": "Caldicellulosiruptor saccharolyticus",
                     "Chlorobium vvibrioforme": "Chlorobium vibrioforme", "Exiguobacterium UNDEF": "Exiguobacterium",
                     "Psychrobacter arcticum": "Psychrobacter arcticus",
                     "Shewanella putefaciens": "Shewanella putrefaciens",
                     "Magnetococcus sp.": "Magnetococcus sp. PL-5-10",
                     "Alkaliphillus metalliredigenes": "Geobacter metallireducens",
                     "Alkalilimnicola ehrlichei": "Alkalilimnicola ehrlichii",
                     "Silicibacter sp.": "Silicibacter sp. 1S17",
                     "Psychrobacter cryopegella": "Psychrobacter cryohalolentis"}

    data = {}
    taxa = set()

    wb = load_workbook(filename = xlsfile, read_only=True)
    for s in wb.get_sheet_names():
        headers = []
        ws=wb[s]
        data[s] = {}
        if verbose:
            if nocolor:
                sys.stderr.write(f"Parsing: {s}\n")
            else:
                sys.stderr.write(f"{bcolors.OKGREEN}Parsing: {s}{bcolors.ENDC}\n")
        for row in ws.rows:
            if not headers:
                for i,j in enumerate(row):
                    if 0 == i and not j.value:
                        headers.append('taxon')
                    if j.value:
                        headers.append(j.value.lower())
                    else:
                        headers.append(None)
                for col in ['gp', 'mr', 'mp', 'blastn']:
                    if col not in headers:
                        if nocolor:
                            sys.stderr.write(f"ERROR: no column named {col} in {s}\n")
                        else:
                            sys.stderr.write(f"{bcolors.WARNING}ERROR:{bcolors.ENDC}: no column named {col} in {s}\n")
                continue
            # ignore empty rows
            if None == row[0].value:
                continue
            # save the first 8 columns
            taxonname = None
            for i,j in enumerate(row):
                if i > 9:
                    break
                if 0 == i:
                    taxonname = j.value
                    spaces = taxonname.count(" ")

                    if spaces > 1:
                        spacesp = taxonname.split(" ")
                        taxonname = spacesp[0] + " " + spacesp[1]

                    if taxonname in correct_names:
                        taxonname = correct_names[taxonname]

                    taxa.add(taxonname)
                    data[s][taxonname] = {}
                else:
                    if not headers[i]:
                        continue
                    if not taxonname:
                        if nocolor:
                            sys.stderr.write(f"FATAL: no taxonomy name\n")
                        else:
                            sys.stderr.write(f"{bcolors.FAIL}FATAL:{bcolors.ENDC}: no taxonomy name\n")
                        sys.exit(-1)
                    data[s][taxonname][headers[i]] = j.value

    return data, taxa

def resolve_taxonomy(tid, conn, verbose=False, nocolor=False):
    """
    Convert the taxonomy id to a tab separated string
    :param tid: the taxonomy object
    :param conn: the database connection
    :param verbose: more output
    :param nocolor: no color ouput
    :return: a string representing the taxonomy
    """

    global taxonomy_str
    if tid in taxonomy_str:
        return taxonomy_str[tid]

    wanted_levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    rnk = ['', '', '', '', '', '', '']
    t, n = get_taxonomy(tid, conn)
    while t.parent != 1 and t.taxid != 1:
        if t.rank in wanted_levels:
            rnk[wanted_levels.index(t.rank)] = n.scientific_name
        t, n = get_taxonomy(t.parent, conn)
    taxonomy_str[tid] = "\t".join(rnk)
    return taxonomy_str[tid]

def add_taxonomy(taxa, verbose=False, nocolor=False):
    """
    Add the taxonomy to our data structure
    :param taxa: the set of taxa to look up
    :param verbose: more output
    :return:
    """

    tids = {}
    taxonomy = {-1 : ['', '', '', '', '', '', '']}

    db = connect_to_db("/data/ncbi/taxonomy.sqlite3")
    if not db:
        sys.stderr.write(f"{bcolors.FAIL}FATAL: could not connect to database{bcolors.ENDC}\n")
        sys.exit(-1)

    for t in taxa:
        tid = get_taxid_for_name(t, db)
        if tid:
            taxonomy[t] = resolve_taxonomy(tid, db)
        else:
            sys.stderr.write(f"{bcolors.WARNING}ERROR: No taxonomy id for {t}{bcolors.ENDC}\n")
            tid = -1
        tids[t] = tid


    return tids, taxonomy

def print_out(data, tids, taxonomy, verbose=False, nocolor=False):
    """
    Print all the data out
    :param data: Out data hash
    :param tids: the taxonomy ids
    :param taxonomy: the taxonomy itself
    :param verbose: more output
    :param nocolor: no color ouput
    :return:
    """

    for s in data:
        for t in data[s]:
            tid = tids[t]
            tax = taxonomy[t]
            tr = []
            for col in ['gp', 'mr', 'mp']:
                tr.append(data[s][t][col])
            truth = "\t".join(map(str, tr))
            sys.stdout.write(f"{s}\t{tid}\t{tax}\t{t}\t{truth}\n")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='excel file', required=True)
    parser.add_argument('-n', help='no color output', action='store_true')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    data, taxa = parse_sheet(args.f, args.v, args.n)
    tids, taxonomy = add_taxonomy(taxa, args.v, args.n)
    print_out(data, tids, taxonomy, args.v, args.n)