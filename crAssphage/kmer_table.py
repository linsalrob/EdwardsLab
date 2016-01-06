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


def select_top(allkmers, percent):
    """
    Select the top percent of kmers

    :param allkmers: A hash of all kmers and their counts
    :type allkmers: hash
    :param percent: The percent to keep
    :type percent: int
    :return: The revised hash
    :rtype: hash
    """


    maxval = max(allkmers.values())
    cutoff = maxval * (1.0 * (100 - percent) / 100)
    allks = {x:allkmers[x] for x in allkmers if allkmers[x] >= cutoff}
    sys.stderr.write("Max value: {}, Cutoff: {}, Length allks: {}\n".format(maxval, cutoff, len(allks)))
    return allks

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="create a kmer table for fecal metagenomes")
    parser.add_argument('-m', help='SRA metadata table from sql query', required=True)
    parser.add_argument('-d', help='directory of jellyfish output files', required=True)
    parser.add_argument('-o', help='output file', required=True)
    parser.add_argument('-p', help='Percent of maximum value to be included. Default = all kmers (-p 100)', default=100, type=int)
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
    allk = {}
    # once we add 100,000 kmers we're only going to add the top n percent (to keep memory usage reasonable!)
    # refactored to store the percentage of kmers as n gets bigger than 2^64!
    addall = True
    for f in os.listdir(args.d):
        if f.endswith('gz'):
            readid = f.split('_')[0]
            if readid in runtype:
                fin = gzip.open(os.path.join(args.d, f), 'rb')
                total = 0
                for l in fin:
                    k, n = l.split()
                    total += int(n)
                fin.close()

                counts[readid] = {}
                if addall:
                    counts[readid][k]=1.0 * int(n)/total
                    allk[k] = allk.get(k, 0) + (1.0 * int(n)/ total)
                else:
                    if k in allk:
                        counts[readid][k] = 1.0 * int(n)/ total
                        allk[k] = allk.get(k, 0) + (1.0 * int(n)/ total)

        if len(allk) > 100000:
            allk = select_top(allk, args.p)
            addall = False


    allks = list(allk.keys())
    allks.sort()

    with open(args.o, 'w') as out:
        for r in counts:
            out.write("{}\t{}".format(r, runtype[r]))
            for k in allks:
                out.write("\t{}".format(counts[r].get(k, 0)))
            out.write("\n")


