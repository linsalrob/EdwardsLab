import os
import re
import sys

from roblib import sequences




def clean_fastq(file1, file2):
    """
    Make a set of cleaned pairs and unpaired reads. If all the reads are paired we do not do anything

    :param file1: first fastq file
    :type file1: str
    :param file2: second fastq file
    :type file2: str
    :return: A string for the output using -1 -2 and -U for unpaired reads
    :rtype: str
    """

    sys.stderr.write("Checking " + file1 + " and " + file2 + "\n")
    seq1 = {}
    seq2 = {}

    for (sid, label, seq, qual) in sequences.stream_fastq(file1):
        sid = re.sub('@', '', sid)
        sid = re.sub('\.[12]$', '', sid)
        sid = re.sub('/[12]$', '', sid)
        seq1[sid] = "@" + label + "\n" + seq + "\n+\n" + qual + "\n"

    for (sid, label, seq, qual) in sequences.stream_fastq(file2):
        sid = re.sub('@', '', sid)
        sid = re.sub('\.[12]$', '', sid)
        sid = re.sub('/[12]$', '', sid)
        seq2[sid] = "@" + label + "\n" + seq + "\n+\n" + qual + "\n"

    seq1set = set(seq1.keys())
    seq2set = set(seq2.keys())

    sys.stderr.write(
        "File 1: " + file1 + " seqs: " + str(len(seq1set)) + " File 2: " + file2 + " seqs: " + str(len(seq2set)) + "\n")

    # are there reads in one but not the other?
    s1unique = seq1set.difference(seq2set)
    s2unique = seq2set.difference(seq1set)

    ret = ' -1 ' + file1 + ' -2 ' + file2

    if len(s1unique) > 0 or len(s2unique) > 0:
        file1 = file1.replace('.gz', '')
        file2 = file2.replace('.gz', '')
        sys.stderr.write("Rewriting " + file1 + " and " + file2 + "\n")
        # we have to make new files


        file1clean = file1.replace('.fastq', '.clean.fastq')
        if file1clean == file1:
            file1clean = file1.replace('.fq', '.clean.fq')
        if file1clean == file1:
            file1clean = file1 + ".clean.fastq"

        file1unique = file1.replace('.fastq', '.unique.fastq')
        if file1unique == file1:
            file1unique = file1.replace('.fq', '.unique.fastq')
        if file1unique == file1:
            file1unique = file1 + '.unique.fastq'

        file2clean = file2.replace('.fastq', '.clean.fastq')
        if file2clean == file2:
            file2clean = file2.replace('.fq', '.clean.fq')
        if file2clean == file2:
            file2clean = file2 + ".clean.fastq"

        file2unique = file2.replace('.fastq', '.unique.fastq')
        if file2unique == file2:
            file2unique = file2.replace('.fq', '.unique.fastq')
        if file2unique == file2:
            file2unique = file2 + '.unique.fastq'

        file1unique = file1unique.replace('.gz', '')
        file2unique = file2unique.replace('.gz', '')
        file1clean = file1clean.replace('.gz', '')
        file2clean = file2clean.replace('.gz', '')


        ret = " -1 " + file1clean + " -2 " + file2clean
        try:
            out1 = open(file1clean, 'w')
            out2 = open(file2clean, 'w')
            for sid in seq1set.intersection(seq2set):
                out1.write(seq1[sid])
                out2.write(seq2[sid])
            out1.close()
            out2.close()
        except IOError as e:
            print "I/O error({0}): {1}".format(e.errno, e.strerror)

        if len(s1unique) > 0:
            ret = ret + " -U " + file1unique
            try:
                out = open(file1unique, 'w')
                for sid in s1unique:
                    out.write(seq1[sid])
                out.close()
            except IOError as e:
                print "I/O error({0}): {1}".format(e.errno, e.strerror)

        if len(s2unique) > 0:
            ret = ret + " -U " + file2unique
            try:
                out = open(file2unique, 'w')
                for sid in s2unique:
                    out.write(seq2[sid])
                out.close()
            except IOError as e:
                print "I/O error({0}): {1}".format(e.errno, e.strerror)
    return ret


import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Given a directory of fastq files, figure out the pairs and check if they are ok. The only requirement is the fastq files have _ between the name and id")
    parser.add_argument('-d', help='Directory of fastq files', required=True)
    args = parser.parse_args()

    files = {}
    for f in os.listdir(args.d):
        sid = f.split('_')[0]
        if sid not in files:
            files[sid] = set()
        files[sid].add(f)

    for s in files:
        if len(files[s]) == 1:
            print(args.d + "/" + files[s].pop())
        elif len(files[s]) == 2:
            outstr = clean_fastq(os.path.join(args.d, files[s].pop()), os.path.join(args.d, files[s].pop()))
            print(outstr)
        else:
            sys.stderr.write("Apparently more than two files for " + s + " ==> " + " ".join(files))
