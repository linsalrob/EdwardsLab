"""
Read a blast file that has gi| numbers and add the taxonomy to the end of the reads. We will add

kingdom phylum class order family genus species
"""

import os
import sys
import time
import argparse
import re
from taxon import taxonomy_hierarchy_as_list, connect_to_db, read_acc_tax_id, acc_to_taxonomy
from roblib import stream_blast_results

default_tax_dir = os.path.join(os.environ['HOME'], "ncbi", "taxonomy", "current")
parser = argparse.ArgumentParser(description="Add taxonomy to a blast output file")
parser.add_argument('-b', help='blast ouput file', required=True)
parser.add_argument('-o', help="Output file (blast results + taxonomy", required=True)
parser.add_argument('-p', help='blast program (blastn, blastp, or blastx (default)', default='blastx')
parser.add_argument('-q', help='Look for taxid of query (default is to use database', action='store_true')
parser.add_argument('-t', help=f'taxonomy directory (default={default_tax_dir})', default=default_tax_dir)
parser.add_argument('-v', help='verbose output', action='store_true')

args = parser.parse_args()

bl_col = 1
if args.q:
    bl_col = 0

dbtype='nucl'
if args.p == 'blastx' or args.p == 'blastp':
    dbtype='prot'

if not os.path.exists(os.path.join(args.t, "taxonomy.sqlite3")):
    print(f"FATAL: taxonomy.sqlite3 not found in {args.t}. Please load the database first", file=sys.stderr)
    sys.exit(1)


conn = connect_to_db(os.path.join(args.t, "taxonomy.sqlite3"))

no_result = ['', '', '', '', '', '', '', '']

results = {}

with open(args.o, 'w') as out:
    for br in stream_blast_results(args.b):
        tid = None
        prot = br.db
        if args.q:
            prot = br.query
        tid, node, name =  acc_to_taxonomy(prot, conn, protein=True, verbose=args.v)
        if not tid:
            if args.v:
                print(f"Warning: No taxid for {prot} found", file=sys.stderr)
            print("\t".join(map(str, br+no_result)), file=out)

        tax = taxonomy_hierarchy_as_list(conn=conn, tid=tid, verbose=args.v)
        print("\t".join(map(str, br+[str(tid)]+tax)), file=out)
