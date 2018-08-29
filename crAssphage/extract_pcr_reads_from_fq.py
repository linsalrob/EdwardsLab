"""
Read the output from extract_pcr_reads -l and pull the left and right reads from a fastq file
"""

import os
import sys
import argparse
import re
from roblib import stream_fastq

def read_readids(regionsf, verbose=True):
    """
    Read the read list from extract_pcr_regions
    :param regionsf: the list
    :param verbose: more output
    :return: a dict of the seq ids and the primer to which they match
    """

    sequences = {}
    with open(regionsf, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            if 3 != len(p):
                continue
            s = re.sub('.\d$', '', p[1])
            if s in sequences and sequences[s] != p[0]:
                sys.stderr.write("ERROR: We had {} for {} and now we have {}\n".format(
                    sequences[s], p[0], s
                ))
            sequences[s] = p[0]
    return sequences

def read_fastqs(fastqdir, fname, seqids, verbose=True):
    """
    Read the fastq files and store the sequences we want to save
    :param fastqdir: the directory with fastq files
    :param fname: the likely filename
    :param seqids: the seqids we want to save
    :param verbose: more output
    :return : a dict of seqids: left, left qual, right, right qual
    """

    seqs = {x:[None, None, None, None] for x in seqids}
    for f in os.listdir(fastqdir):
        if fname in f:
            for seqid, header, seq, qualscores in stream_fastq(os.path.join(fastqdir, f)):
                s = re.sub('.\d$', '', seqid)
                if s in seqids:
                    if seqid.endswith('.1'):
                        seqs[s][0] = seq
                        seqs[s][1] = qualscores
                    elif seqid.endswith('.2'):
                        seqs[s][2] = seq
                        seqs[s][3] = qualscores
                    else:
                        sys.stderr.write("ERR: Not sure about sequence ID {}\n".format(seqid))
    return seqs

def write_fastq(seqids, seqs, fname, verbose=False):
    """
    Write the fastq files
    :param seqids: sequence ids and the primers to which they belong
    :param seqs: sequences
    :param fname: the likely filename
    :param verbose: more output
    :return: nothing
    """
    for primer in ['A', 'B', 'C']:
        lft = open("{}_Primer{}_1.fastq".format(fname, primer), 'w')
        rht = open("{}_Primer{}_2.fastq".format(fname, primer), 'w')
        sng = open("{}_Primer{}_single.fastq".format(fname, primer), 'w')
        for s in seqids:
            if seqids[s] == 'Primer{}'.format(primer):
                if seqs[s][0] and seqs[s][2]:
                    lft.write("@{}.1\n{}\n+\n{}\n".format(s, seqs[s][0], seqs[s][1]))
                    rht.write("@{}.2\n{}\n+\n{}\n".format(s, seqs[s][2], seqs[s][3]))
                elif seqs[s][0]:
                    sng.write("@{}\n{}\n+\n{}\n".format(s, seqs[s][0], seqs[s][1]))
                elif seqs[s][2]:
                    sng.write("@{}\n{}\n+\n{}\n".format(s, seqs[s][2], seqs[s][3]))
                else:
                    sys.stderr.write("No sequences for {}\n".format(s))
        lft.close()
        rht.close()
        sng.close()










if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="make L and R fastq files from extract_pcr_regions output")
    parser.add_argument('-l', help='output from extract_pcr_regions with the -l flag', required=True)
    parser.add_argument('-d', help='directory with fastq files of the data. We assume that the file name is the same for -l and the fastq file except for the last part after the .', required=True)
    parser.add_argument('-o', help='output directory', required=True)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    basefile = args.l
    if '/' in basefile:
        basefile = basefile.split('/')[-1]
    basefile = basefile.split('.')[0]

    sys.stderr.write("Assuming your base file name is {}\n".format(basefile))

    seqids = read_readids(args.l, args.v)
    seqs = read_fastqs(args.d, basefile, seqids, args.v)
    write_fastq(seqids, seqs, basefile, args.v)