"""
Read a fastq file and a time file and split based on the time time stamp

Our time stamp file should have the following information:

sample ID\tStart Time

We assume that all times are contiguous, so if you want to include a gap you'll need to add that with a name

Time should be in one of the formats supported by `dateutil.parser <http://dateutil.readthedocs.io/en/stable/parser.html>`

e.g.:
 - 2018-06-13T09:12:57Z
 -  2018-06-13T19:12:57AEST

e.g.:
F6      2018-06-13T12:23:00AEST
B10     2018-06-13T14:41:00AEST
D8      2018-06-13T16:41:00AEST
C7      2018-06-13T18:42:00AEST
B12     2018-06-13T20:40:00AEST
G1      2018-06-14T08:56:00AEST
D10     2018-06-14T10:55:00AEST



"""

import os
import sys
import argparse
import re

from dateutil.parser import parse
from datetime import timedelta

from fastq import stream_fastq


def write_fastq(fqf, outs, outdir):
    """
    Write the sequences to a set of fastq files
    :param fqf: the input fastq file with the original sequences
    :param outs: the sets of sequences for each id
    :param outdir: the output directory to write the sequences to
    :return:
    """

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    outputs = {}
    for o in outs:
        outputs[o] = open(os.path.join(outdir, o + ".fastq"), 'w')

    remapped = {}
    for o in outs:
        for seqid in outs[o]:
            remapped[seqid] = o

    for seqid, header, seq, qualscores in stream_fastq(fqf):
        if seqid not in remapped:
            sys.stderr.write("Error: found sequence {} that we don't know where to write{}\n".format(seqid))
        outputs[remapped[seqid]].write("@{}\n{}\n+\n{}\n".format(header, seq, qualscores))

    for o in outputs:
        outputs[o].close()


def split_fastq(fqf, times):
    """
    Split the fastq file based on the times and dates
    :param fqf: fastq file to parse
    :param times: dictionary of times
    :return: a dictionary of time ids and the list of sequences in that time
    """

    seqs = {"other" : set(), "before" : set(), "after" : set()}
    alltimes = []

    for t in times:
        seqs[t] = set()
        alltimes.append(times[t][0])
        alltimes.append(times[t][1])

    earliest = min(alltimes)
    latest   = max(alltimes)

    newest = None
    newestseq = None
    oldest = None
    oldestseq = None

    for seqid, header, seq, qualscores in stream_fastq(fqf):
        m = re.search("start_time=([\w\:\-]+)", header)
        if not m:
            sys.stderr.write("No start time was detected in {}\n".format(header))
            continue
        try:
            seqstart = parse(m.groups()[0])
        except ValueError as v:
            sys.stderr.write("Can't parse date time from: {}\n".format(m.groups()[0]))
            continue

        if seqstart < earliest:
            seqs['before'].add(seqid)
            continue
        if seqstart > latest:
            seqs['after'].add(seqid)
            continue

        if not newest or seqstart < newest:
            newest = seqstart
            newestseq = seqid

        if not oldest or seqstart > oldest:
            oldest = seqstart
            oldestseq = seqid

        added = False
        for t in times:
            if seqstart > times[t][0] and seqstart <= times[t][1]:
                added = True
                seqs[t].add(seqid)
                break
        if not added:
            seqs['other'].add(seqid)

    sys.stderr.write("Newest sequence: {} at {}\nOldest sequence: {} at {}\n".format(
        newestseq, newest, oldestseq, oldest
    ))

    return seqs

def parse_times(timefile, ztoffset):
    """
    Parse the times from the time separation file
    :param timefile: the file to parse
    :param ztoffset: the difference from zulu time
    :return: a dict of IDs and times
    """

    times = {}
    lastid = None
    with open(timefile, 'r') as fin:
        for l in fin:
            p=l.strip().split("\t")
            try:
                starttime = parse(p[1])
            except:
                sys.stderr.write("Error: could not parse start time from {}\n".format(p[1]))
                continue
            if ztoffset:
                startime = starttime + timedelta(hours=ztoffset)

            times[p[0]] = [starttime, 0]
            if lastid:
                times[lastid][1] = starttime
            lastid = p[0]

    times[lastid][1] = times[lastid][0] + timedelta(hours=48)

    for t in times:
        sys.stderr.write("Time: {} From: {} To: {}\n".format(t, times[t][0], times[t][1]))

    return times



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Deconvolute a fastq file")
    parser.add_argument('-f', help='fastq file to read', required=True)
    parser.add_argument('-t', help='timestamp file', required=True)
    parser.add_argument('-z', help='offset from zulu time. This number will be added to Z so use -8 on the west coast of US. Default: use Z time', type=int, default=0)
    parser.add_argument('-o', help='output directory to write fastq files to. If not provided sequences not written')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    times = parse_times(args.t, args.z)
    seqs  = split_fastq(args.f, times)

    if args.o:
        write_fastq(args.f, seqs, args.o)

    for t in seqs:
        sys.stdout.write("{}\t{}\n".format(t, "; ".join(seqs[t])))


