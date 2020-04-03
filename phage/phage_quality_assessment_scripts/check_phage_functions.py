
"""
Test the functions of the phages.
"""

import os
import sys
import argparse
from pppf_databases import connect_to_db, disconnect
from pppf_clusters import proteinid_to_function
from roblib import is_hypothetical

def check_phage_functions(sample, blastfile, outputfile, clusterdb):
    """
    Count how many proteins have hypothetical functions
    """

    out = open(outputfile, 'w')
    if not os.path.exists(clusterdb):
        out.close()
        return

    phage_cluster_db = connect_to_db(clusterdb)
    hypo = 0
    nonhypo = 0
    with open(blastfile, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            fn = proteinid_to_function(p[1], phage_cluster_db)
            if is_hypothetical(fn):
                hypo += 1
            else:
                nonhypo += 1
    out.write(f"{sample}\tHypothetical proteins\t")
    out.write("[Hypothetical, Non-hypothetical, Fraction hypothetical]\t")
    out.write(f"{hypo}\t{nonhypo}\t{hypo / (hypo + nonhypo)}\n")
    out.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-s', help='sample name used in output', required=True)
    parser.add_argument('-b', help='blast m8 file', required=True)
    parser.add_argument('-o', help='output file', required=True)
    parser.add_argument('-c', help='phage clusters dv', required=True)
    args = parser.parse_args()

    check_phage_functions(args.s, args.b, args.o, args.c)