"""
Yet another attempt to separate sequences for prophages, This uses the coordinates
output from phispy 20190920. Note that there is an additional tab here,
that we currently ignore!
"""

import os
import sys
import argparse
from operator import itemgetter
from roblib import read_fasta, bcolors

def read_phage_locations(locf, verbose=False):
    """
    Read the phage locations. Expects tab separated text. The format is 
    [pp#, blank, contig, start, stop] and then potentially repeat information
    [r1start, r1end, r2start, r2end, sequence, information]

    :param locf: the locations file
    :param verbose: more output
    :return: a dict of all the contigs that have prophages, and then tuples of [start, stop]
    """

    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}READING {locf}{bcolors.ENDC}\n")
    contigs = {}
    with open(locf, 'r') as f:
        for l in f:
            p=l.strip().split("\t")
            try:
                c = p[2] # contig
                s = int(p[3]) # start
                e = int(p[4]) # end
            except:
                sys.stderr.write(f"{bcolors.FATAL}ERROR PARSING {l} in {locf}\n{bcolors.ENDC}")
                sys.exit(-1)

            if s > e:
                (e, s) = (s, e)
            if c not in contigs:
                contigs[c] = []
            contigs[c].append((s, e))
    return contigs


def read_sequence(conf, verbose=False):
    """
    Read the contigs file for this genome and return it
    :param conf: the contigs file
    :param verbose:
    :return: a dict of contig/sequence
    """
    
    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}READING {conf}{bcolors.ENDC}\n")
    return read_fasta(conf, whole_id = False)

def write_seqs(seqs, locations, phagef, nonphagef, minlen, verbose=False):
    """
    Write the sequences out
    :param phagef: file to write the prophages to
    :param nonphagef: file to write the nonprophages to
    :param seqs: DNA sequences
    :param locations: locations of phage regions
    :param minlen: minimum length of sequence to write out
    :param verbose: more output
    :return:
    """
    
    if verbose:
        sys.stderr.write(f"{bcolors.GREEN}Writing sequences to {phagef} and {nonphagef}{bcolors.ENDC}\n")

    with open(phagef, 'w') as phageout:
        with open(nonphagef, 'w') as nonphageout:
            for s in seqs:
                if s in locations:
                    ses = sorted(locations[s], key=itemgetter(0))
                    posn = 0
                    for start, end in ses:
                        if verbose:
                            sys.stderr.write(f"{bcolors.BLUE}{posn} to {start-1} and then {start} - {end} ")
                            sys.stderr.write(f"seq len: {len(seqs[s])}{bcolors.ENDC}\n")

                        seqid = '_'.join(map(str, [s, posn, start-1]))
                        seqstr = seqs[s][posn:start]
                        if len(seqstr) > minlen:
                            nonphageout.write(f">{seqid}\n{seqstr}\n")

                        seqid = '_'.join(map(str, [s, start, end]))
                        seqstr = seqs[s][start:end]
                        if len(seqstr) > minlen:
                            phageout.write(f">{seqid}\n{seqstr}\n")
                        posn=end+1
                    
                    # write the rest of the sequence
                    seqid = '_'.join(map(str, [s, posn, len(seqs[s])]))
                    seqstr = seqs[s][posn:]
                    if len(seqstr) > minlen:
                        nonphageout.write(f">{seqid}\n{seqstr}\n")
                else:
                    if verbose:
                        sys.stderr.write(f"{bcolors.PINK}|{s} not in locations{bcolors.ENDC}\n")
                    nonphageout.write(">{}\n{}\n".format("_".join(map(str, [s, 0, len(seqs[s])])),
                                                       seqs[s]))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Separate prophage and non prophage DNA")
    parser.add_argument('-l', help='locations in tsv format file', required=True)
    parser.add_argument('-c', help='contigs file', required=True)
    parser.add_argument('-m', help='minimum length of sequence to save at the end of a contig, (default=100bp)', default=100, type=int)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    for d in ['prophage_seqs', 'nonprophage_seqs']:
        if not os.path.exists(d):
            os.mkdir(d)

    '''
    This may be useful for converting to a whole directory
    
    for f in os.listdir(args.l):
        if os.path.exists(os.path.join(args.c, f"{f}.fna")):
            loc = read_phage_locations(os.path.join(args.l, f))
            dna = read_sequence(os.path.join(args.c, f"{f}.fna"))
            write_seqs(dna, loc, os.path.join("prophage_seqs", f), os.path.join("nonprophage_seqs", f))
    '''

    loc = read_phage_locations(args.l, args.v)
    dna = read_sequence(args.c, args.v)
    outputfile = args.c.split("/")[-1]
    if '.fasta' in outputfile:
        bout = outputfile.replace('.fasta', '.bacteria.fasta') 
        pout = outputfile.replace('.fasta', '.prophage.fasta')
    elif '.fna' in outputfile:
        bout = outputfile.replace('.fna', '.bacteria.fna') 
        pout = outputfile.replace('.fna', '.prophage.fna')
    else:
        if verbose:
            sys.stderr.write(f"{bcolors.YELLOW}WARNING: {outputfile} does not contain fasta or fna{bcolors.ENDC}\n")
        bout = outputfile + ".bacteria"
        pout = outputfile + ".prophages"

    write_seqs(dna, loc, os.path.join("prophage_seqs", pout), os.path.join("nonprophage_seqs", bout), args.m, args.v) 
