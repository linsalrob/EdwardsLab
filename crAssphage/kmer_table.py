"""
Make a data table by combining the fecal run ids extracted from the SQL
with the 21-mers counted by jellyfish.

We make the run ids table by using this sql query (it takes a while to run!)

sqlite3 SRAmetadb.sqlite "select run_accession, e.experiment_accession, e.library_strategy, s.study_type from run
as r left join experiment as e on r.experiment_accession = e.experiment_accession left join study as s on
s.study_accession = e.study_accession where e.experiment_accession in (select e.experiment_accession from experiment
as e where e.study_accession in (select study_accession from study where study_type like 'metagenomics'));"
> metagenomics_info.txt

We count the jellyfish 21-mers in the usual way. Now we just need to munge the two tables!

"""
import gzip
import os
import sys
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="create a kmer table for fecal metagenomes")
    parser.add_argument('-m', help='SRA metadata table from sql query', required=True)
    parser.add_argument('-d', help='directory of jellyfish output files', required=True)
    parser.add_argument('-o', help='output file', required=True)
    args = parser.parse_args()

    runtype = {}
    with open(args.m, 'r') as min:
        for l in min:
            if '\t' in l:
                p = l.split('\t')
            elif '|' in l:
                p = l.split('|')
            else:
                sys.exit("Ooops. What is the field delimiter in {}\n".format(args.m))
            runtype[p[0]]=p[2]

    counts = {}
    allk = set()
    for f in os.listdir(args.d):
        if f.endswith('gz'):
            readid = f.split('_')[0]
            if readid in runtype:
                counts[readid] = {}
                fin = gzip.open(os.path.join(args.d, f), 'rb')
                for l in fin:
                    k, n = l.split()
                    counts[readid][k]=n
                    allk.add(k)

    allks = list(allk)
    allks.sort()
    with open(args.o, 'w') as out:
        for r in counts:
            sys.stderr.write("{}\n".format(r))
            out.write("{}\t{}".format(r, runtype[r]))
            for k in allks:
                out.write("\t{}".format(counts[r].get(k, 0)))
            out.write("\n")


