"""
Check a fasta file to make sure it is correct
"""

import os
import sys
import argparse



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Check the integrity of a fasta file")
    parser.add_argument('-f', help='fasta file to check', required=True)
    parser.add_argument('-o', help='output file to write', required=True)
    parser.add_argument('-s', help='strip new lines from sequence (default: keep new lines)', action="store_true")
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()
    
    seqid = None
    seq = ""

    with open(args.f, 'r') as fin:
        with open(args.o, 'w') as out:
            for l in fin:
                l = l.strip()
                if not l:
                    continue

                if l.startswith(">"):
                    if seqid:
                        if seq:
                            out.write(seqid)
                            out.write(seq)
                        elif args.v:
                            sys.stderr.write("There is no sequence for {}. Skipped\n".format(seqid))
                    seqid = l + "\n"
                    seq = ""
                    continue

                if ">" in l:
                    # line contains some sequence and then a new header line
                    tmpseq = l[:l.index('>')]
                    tmpid = l[l.index('>'):]
                    seq += tmpseq + "\n"
                    out.write(seqid)
                    out.write(seq)
                    seqid = tmpid + "\n"
                    seq = ""
                    continue

                if args.s:
                    seq += l
                else:
                    seq += l + "\n"

            if seqid:
                if seq:
                    out.write(seqid)
                    out.write(seq)
                elif args.v:
                    sys.stderr.write("There is no sequence for {}. Skipped\n".format(seqid))