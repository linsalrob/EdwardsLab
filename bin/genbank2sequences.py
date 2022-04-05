#!/usr/bin/env python
"""
Convert a genbank file to sequences
"""

import os
import sys
import gzip
import argparse
from roblib import genbank_to_faa, genbank_to_fna, genbank_to_orfs, genbank_to_ptt, genbank_to_functions, genbank_seqio
from roblib import genbank
from Bio import SeqIO

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-g', '--genbank', help='genbank file', required=True)
    parser.add_argument('-c', '--complex', help='complex identifier line', action='store_true')
    parser.add_argument('-a', '--aminoacids', help="output file for the amino acid sequences (.faa will be appended)")
    parser.add_argument('-n', '--nucleotide', help='output file for nucleotide sequence (.fna will be appended)')
    parser.add_argument('-p', '--ptt', help='output file for the ptt protein table')
    parser.add_argument('-o', '--orfs', help='output file for orfs (.orfs will be appended)')
    parser.add_argument('-f', '--functions', help='output file for two column table of [protein id, function]')
    parser.add_argument('-i', '--seqid', help='Only output these sequence ID(s) [multiple -i allowed]',
                        action='append')
    parser.add_argument('--gff3', help="Output gff3 format (experimental)")
    parser.add_argument('--phage_finder', help='make a phage finder file')
    parser.add_argument('--separate',  action='store_true',
                        help='separate output into different files (with no other options just output gbk).')
    parser.add_argument('-z', '--zip', help='gzip compress the output. Experimental and may not work with everything!',
                        action='store_true')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not os.path.exists(args.genbank):
        sys.stderr.write(f"FATAL: {args.genbank} does not exist. Please check the file path and try again!")
        sys.exit(1)

    if args.seqid and not args.separate:
        sys.stderr.write("-i was provided, so requiring to separate files (--separate assumed)\n")
        args.separate = True

    did = False
    if args.nucleotide:
        if args.separate:
            lastid = None
            out = None
            for sid, seq in genbank_to_fna(args.genbank, args.complex):
                if args.seqid and sid not in args.seqid:
                    if args.v:
                        sys.stderr.write(f"Skipped {sid} not provided in -i options\n")
                    continue
                if sid != lastid:
                    if out:
                        out.close()
                    out = open(f"{args.nucleotide}.{sid}.fna", 'w')
                    lastid = sid
                out.write(f">{sid}\n{seq}\n")
            if out:
                out.close()
        else:
            with open(f"{args.nucleotide}.fna", 'w') as out:
                for sid, seq in genbank_to_fna(args.genbank, args.complex):
                    out.write(f">{sid}\n{seq}\n")
        did = True

    if args.aminoacids:
        if args.separate:
            lastid = None
            out = None
            for seqid, sid, seq in genbank_to_faa(args.genbank, args.complex, args.v):
                if args.seqid and sid not in args.seqid:
                    if args.v:
                        sys.stderr.write(f"Skipped {seqid} not provided in -i options\n")
                    continue
                if seqid != lastid:
                    if out:
                        out.close()
                    out = open(f"{args.aminoacids}.{seqid}.faa", 'w')
                    lastid = seqid
                out.write(f">{sid}\n{seq}\n")
            if out:
                out.close()
        else:
            with open(f"{args.aminoacids}.faa", 'w') as out:
                for seqid, sid, seq in genbank_to_faa(args.genbank, args.complex, args.v):
                    out.write(f">{sid}\n{seq}\n")
        did = True

    if args.orfs:
        if args.separate:
            lastid = None
            out = None
            for seqid, sid, seq in genbank_to_orfs(args.genbank, args.complex, args.v):
                if args.seqid and sid not in args.seqid:
                    if args.v:
                        sys.stderr.write(f"Skipped {seqid} not provided in -i options\n")
                    continue
                if seqid != lastid:
                    if out:
                        out.close()
                    out = open(f"{args.orfs}.{seqid}.orfs", 'w')
                    lastid = seqid
                out.write(f">{sid}\n{seq}\n")
            if out:
                out.close()
        else:
            with open(f"{args.orfs}.orfs", 'w') as out:
                for seqid, sid, seq in genbank_to_orfs(args.genbank, args.complex, args.v):
                    out.write(f">{sid}\n{seq}\n")
        did = True

    if args.ptt:
        r = genbank_to_ptt(args.genbank, False, args.v)
        with open(args.ptt, 'w') as out:
            for l in r:
                out.write("\t".join(map(str, l)))
                out.write("\n")
        did = True

    if args.functions:
        try:
            if args.zip:
                out = gzip.open(f"{args.functions}.gz", 'wt')
            else:
                out = open(args.functions, 'w')
            for sid, pid, prod in genbank_to_functions(args.genbank, True, args.v):
                out.write(f"{sid}\t{pid}\t{prod}\n")
            did = True
            out.close()
        except IOError as e:
            sys.stderr.write(f"There was an error writing to {args.functions}: {e}\n")
            sys.exit(1)

    if args.phage_finder:
        with open(args.phage_finder, 'w') as out:
            for tple in genbank.genbank_to_phage_finder(args.genbank, args.v):
                out.write("\t".join(map(str, tple)) + "\n")
        did = True

    if args.gff3:
        genbank_to_gff(args.genbank, args.gff3, args.v)
        did = True

    if not did and args.separate:
        lastid = None
        out = None
        for seq in genbank_seqio(args.genbank):
            if args.seqid and seq.id not in args.seqid:
                if args.v:
                    sys.stderr.write(f"Skipped {seq.id} not provided in -i options\n")
                continue
            out = open(f"{seq.id}.gbk", 'w')
            SeqIO.write(seq, out, 'genbank')
            out.close()
        did = True
    
    if not did:
        sys.stderr.write("Please provide either a -n, -a, -o, -p, -f, --gff3 output file! (or all)\n")
