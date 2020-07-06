"""
Rub hmmer and print the results.

This is really test code for PhiSpy, but may be useful elsewhere!

Now it has turned into timing tests!
"""

import os
import sys
import argparse
from roblib import colours, genbank_seqio, feature_id
import subprocess
from Bio import SearchIO
from io import StringIO
import timeit
from tempfile import NamedTemporaryFile

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


def run_hmmscan_oat(gbkf, hmmf):

    allhits = {}
    hitcount = 0
    rescount = 0
    for seq in genbank_seqio(gbkf):
        for feat in seq.features:
            if feat.type != 'CDS':
                continue
            aa = ""
            if 'translation' in feat.qualifiers:
                aa = feat.qualifiers['translation'][0]
            else:
                aa = str(feat.extract(seq).translate().seq)

            search = subprocess.Popen(["hmmscan", '--cpu', '6', '-E', '1e-10', '--domE', '1e-5', '--noali', hmmf, '-'],
                                      stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            hmmresult = search.communicate(input=f">{feat.id}\n{aa}\n".encode())[0]

            results = SearchIO.parse(StringIO(hmmresult.decode()), 'hmmer3-text')
            for res in results:
                allhits[res.id] = {}
                rescount += 1
                for hit in res:
                    allhits[res.id][hit.id] = hit.evalue
                    hitcount += 1

    print(f"Using hmmscan and streaming one at time there were {rescount} results and {hitcount} hits, and our dict has {len(allhits)} entries")


def run_hmmscan_aao(gbkf, hmmf):


    allhits = {}
    hitcount = 0
    rescount = 0
    for seq in genbank_seqio(gbkf):
        prots = []
        for feat in seq.features:
            if feat.type != 'CDS':
                continue
            aa = ""
            if 'translation' in feat.qualifiers:
                aa = feat.qualifiers['translation'][0]
            else:
                aa = str(feat.extract(seq).translate().seq)
            myid = feature_id(seq, feat)
            prots.append(f">{myid}\n{aa}")

        search = subprocess.Popen(["hmmscan", '--cpu', '6', '-E', '1e-10', '--domE', '1e-5', '--noali', hmmf, '-'],
                                  stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        hmmresult = search.communicate(input="\n".join(prots).encode())[0]

        results = SearchIO.parse(StringIO(hmmresult.decode()), 'hmmer3-text')
        for res in results:
            allhits[res.id] = {}
            rescount += 1
            for hit in res:
                allhits[res.id][hit.id] = hit.evalue
                hitcount += 1

    print(f"Using hmmscan and streaming all at once there were {rescount} results and {hitcount} hits, and our dict has {len(allhits)} entries")

def hmmscan_print_then_run(gbkf, hmmf):
    aaout = NamedTemporaryFile(mode='w+t', delete=False)
    sys.stderr.write(f"Writing the amino acids to {aaout.name}\n")
    aaout.seek(0)
    for seq in genbank_seqio(gbkf):
        for feat in seq.features:
            if feat.type != 'CDS':
                continue
            aa = ""
            if 'translation' in feat.qualifiers:
                aa = feat.qualifiers['translation'][0]
            else:
                aa = str(feat.extract(seq).translate().seq)
            myid = feature_id(seq, feat)
            aaout.write(f">{myid}\n{aa}\n")

    aaout.close()

    hmmout = NamedTemporaryFile(mode='w+t', delete=False)
    hmmout.close()
    sys.stderr.write(f"hmmscan output will be in {hmmout.name}\n")
    try:
        search = subprocess.run(["hmmscan", '--cpu', '6', '-E', '1e-10', '--domE', '1e-5', '--noali', '-o', hmmout.name, hmmf, aaout.name])
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Error running hmmscan:\n{e}\n")
        sys.exit(-1)


    results = SearchIO.parse(hmmout.name, 'hmmer3-text')
    allhits = {}
    hitcount = 0
    rescount = 0
    for res in results:
        allhits[res.id] = {}
        rescount += 1
        for hit in res:
            allhits[res.id][hit.id] = hit.evalue
            hitcount += 1

    print(f"Using hmmscan and tempfiles there were {rescount} results and {hitcount} hits, and our dict has {len(allhits)} entries")


def hmmsearch_print_then_run(gbkf, hmmf):
    aaout = NamedTemporaryFile(mode='w+t', delete=False)
    sys.stderr.write(f"Writing the amino acids to {aaout.name}\n")
    aaout.seek(0)
    for seq in genbank_seqio(gbkf):
        for feat in seq.features:
            if feat.type != 'CDS':
                continue
            aa = ""
            if 'translation' in feat.qualifiers:
                aa = feat.qualifiers['translation'][0]
            else:
                aa = str(feat.extract(seq).translate().seq)
            myid = feature_id(seq, feat)
            aaout.write(f">{myid}\n{aa}\n")

    aaout.close()

    hmmout = NamedTemporaryFile(mode='w+t', delete=False)
    hmmout.close()
    sys.stderr.write(f"hmmsearch output will be in {hmmout.name}\n")
    try:
        search = subprocess.run(["hmmsearch", '--cpu', '6', '-E', '1e-10', '--domE', '1e-5', '--noali', '-o', hmmout.name, hmmf, aaout.name])
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Error running hmmscan:\n{e}\n")
        sys.exit(-1)


    results = SearchIO.parse(hmmout.name, 'hmmer3-text')
    allhits = {}
    hitcount = 0
    rescount = 0
    for res in results:
        allhits[res.id] = {}
        rescount += 1
        for hit in res:
            allhits[res.id][hit.id] = hit.evalue
            hitcount += 1

    print(f"Using hmmsearch and tempfiles there were {rescount} results and {hitcount} hits, and our dict has {len(allhits)} entries")



def hmmsearch_print_then_parse(gbkf, hmmf):
    aaout = NamedTemporaryFile(mode='w+t', delete=False)
    sys.stderr.write(f"Writing the amino acids to {aaout.name}\n")
    aaout.seek(0)
    for seq in genbank_seqio(gbkf):
        for feat in seq.features:
            if feat.type != 'CDS':
                continue
            aa = ""
            if 'translation' in feat.qualifiers:
                aa = feat.qualifiers['translation'][0]
            else:
                aa = str(feat.extract(seq).translate().seq)
            myid = feature_id(seq, feat)
            aaout.write(f">{myid}\n{aa}\n")

    aaout.close()

    sys.stderr.write("Searching\n")
    try:
        search = subprocess.Popen(["hmmsearch", '--cpu', '6', '-E', '1e-10', '--domE', '1e-5', '--noali', hmmf, aaout.name], stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Error running hmmscan:\n{e}\n")
        sys.exit(-1)

    sys.stderr.write("Parsing\n")
    hmmresult = search.communicate()[0]
    results = SearchIO.parse(StringIO(hmmresult.decode()), 'hmmer3-text')
    allhits = {}
    hitcount = 0
    rescount = 0
    for res in results:
        allhits[res.id] = {}

        rescount += 1
        for hit in res:
            allhits[res.id][hit.id] = hit.evalue
            # print(f"Result: {res.id}: Hit: {hit.id} Eval: {hit.evalue}")
            hitcount += 1

    print(f"Using hmmsearch and tempfiles there were {rescount} results and {hitcount} hits, and our dict has {len(allhits)} entries")



def stream_hmmsearch(gbkf, hmmf):
    """
    NOTE: THIS DOES NOT WORK!!

    You can't stream against hmmsearch as it can't rewind the sequences. You either need to use hmmscan (slow) or a temp file
    (fast!)

    :param gbkf:
    :param hmmf:
    :return:
    """
    allhits = {}
    hitcount = 0
    rescount = 0
    for seq in genbank_seqio(gbkf):
        prots = []
        for feat in seq.features:
            if feat.type != 'CDS':
                continue
            aa = ""
            if 'translation' in feat.qualifiers:
                aa = feat.qualifiers['translation'][0]
            else:
                aa = str(feat.extract(seq).translate().seq)
            myid = feature_id(seq, feat)
            prots.append(f">{myid}\n{aa}")


        search = subprocess.Popen(["hmmsearch", '--cpu', '6', '-E', '1e-10', '--domE', '1e-5', '--noali', hmmf, '-'],
                                  stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        hmmresult = search.communicate(input="\n".join(prots).encode())[0]

        results = SearchIO.parse(StringIO(hmmresult.decode()), 'hmmer3-text')
        for res in results:
            allhits[res.id] = {}
            rescount += 1
            for hit in res:
                allhits[res.id][hit.id] = hit.evalue
                hitcount += 1

    print(
        f"Using hmmsearch and streaming all at once there were {rescount} results and {hitcount} hits, and our dict has {len(allhits)} entries")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-g', help='genbank file', required=True)
    parser.add_argument('-m', help='hmmer file', required=True)
    args = parser.parse_args()


    print("Timing Search then parse")
    print(timeit.timeit('hmmsearch_print_then_parse(args.g, args.m)', setup="from __main__ import hmmsearch_print_then_parse, args", number=3))
    # takes about 600 seconds

    print("Timing hmmsearch using temporary files")
    print(timeit.timeit('hmmsearch_print_then_run(args.g, args.m)', setup="from __main__ import hmmsearch_print_then_run, args", number=3))
    # takes about 600 seconds
    exit(0)

    print("Timing hmmsearch")
    print(timeit.timeit('stream_hmmsearch(args.g, args.m)', setup="from __main__ import stream_hmmsearch, args", number=3))
    # Does not work (see note above)

    print("Timing run all at once")
    print(timeit.timeit('run_hmmscan_aao(args.g, args.m)', setup="from __main__ import run_hmmscan_aao, args", number=3))
    # takes about 3000 seconds

    print("Timing run one at a time")
    print(timeit.timeit('run_hmmscan_oat(args.g, args.m)', setup="from __main__ import run_hmmscan_oat, args", number=3))
    # takes about 3000 seconds

    print("Timing hmmscan using temporary files")
    print(timeit.timeit('hmmscan_print_then_run(args.g, args.m)', setup="from __main__ import hmmscan_print_then_run, args", number=3))
    # takes about 3000 seconds

    print("Timing hmmsearch using temporary files")
    print(timeit.timeit('hmmsearch_print_then_run(args.g, args.m)', setup="from __main__ import hmmsearch_print_then_run, args", number=3))
    # takes about 600 seconds

