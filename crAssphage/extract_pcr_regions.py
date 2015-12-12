import argparse
import os
import sys
import pysam

__author__ = 'Rob Edwards'

"""
Extract the PCR regions from indexed bam alignments. These are regions that we have asked people to pull out, and
we will extract the sequences from those regions

Current regions:


Primer sequences:			Position on JQ995537 (the original genbank record)
Primer A:
	Fwd: CTGATAGTATGATTGGTAAT	Position 25634 .. 25653
	Rev: ATAAGTTCTCCAACTATCTT	Position complement(26945 .. 26964)

Primer B:
	Fwd: CCAGTATCTCCATAAGCATC	Position 33709 .. 33728
	Rev: GTGAGGGCGGAATAGCTA		Position complement(35045 .. 35062)

Primer C:
	Fwd: GCAACAGGAGTAGTAAAATCTC	Position 43820 .. 43841
	Rev: GCTCCTGTTAATCCTGATGTTA	Position complement(45036 .. 45057)

"""

locations = {
    'JQ995537': {
        'PrimerA': (25634, 26964),
        'PrimerB': (33709, 35062),
        'PrimerC': (43820, 45057)
    }
}


def print_query_regions(bam):
    """
    Print the regions from the query sequences that overlap our regions of interest

    :param bam: The bam object from pysam
    :type bam: pysam.AlignmentFile
    :return:
    :rtype:
    """

    for template in locations:
        for primer in locations[template]:
            start, end = locations[template][primer]
            for read in bam.fetch(reference=template, start=start, end=end):
                # this is an AlignedSegment: http://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedSegment
                # sys.stderr.write("Primer: {} ({} .. {}). Found a region for {} ({} .. {}) -> ({} .. {})\n".format(
                #     primer, start, end, read.query_name, read.query_alignment_start, read.query_alignment_end,
                #     read.reference_start, read.reference_end
                # ))

                # this checks for sequences that overlap the start and end (none do in the Ondrej data set
                # if read.reference_start <= start and read.reference_end >= stop:
                #     sys.stderr.write("Primer: {} ({} .. {}). Found a region for {} ({} .. {}) -> ({} .. {})\n".format(
                #         primer, start, stop, read.query_name, read.query_alignment_start, read.query_alignment_end,
                #         read.reference_start, read.reference_end
                #     ))

                # get just the sequence that maps to the region
                seq = read.query_sequence
                beg_offset = None
                end_offset = None
                if read.reference_start < start:
                    beg_offset = start - read.reference_start - 1
                if read.reference_end > end:
                    end_offset = len(seq) - (read.reference_end - end)

                if beg_offset and end_offset:
                    seq = seq[beg_offset:end_offset]
                elif beg_offset:
                    seq = seq[beg_offset:]
                elif end_offset:
                    seq = seq[:end_offset]

                print(">{} {} {} {}\n{}".format(read.query_name, primer, read.reference_start, read.reference_end, seq))

def print_pileup(bam):
    """
    Print the information about the pileup

    :param bam: the bam object from pysam
    :type bam: pysam.AlignmentFile
    :return:
    :rtype:
    """

    for template in locations:
        for primer in locations[template]:
            start, end  = locations[template][primer]
            for p in bam.pileup(reference=template, start=start, end=end, truncate=True):
                bases = {}
                for pilups in p.pileups:
                    if pilups.query_position:
                        bp = pilups.alignment.query_sequence[pilups.query_position]
                    else:
                        bp = '-'
                    bases[bp] = bases.get(bp, 0) + 1
                sys.stdout.write("{} : {} -> {}\n".format(p.reference_name, p.reference_pos, str(bases)))

def print_consensus(bam):
    """
    Print the consensus sequence about the pileup

    :param bam: the bam object from pysam
    :type bam: pysam.AlignmentFile
    :return:
    :rtype:
    """
    for template in locations:
        for primer in locations[template]:
            start, end = locations[template][primer]
            cseq = []
            for p in bam.pileup(reference=template, start=start, end=end, truncate=True):
                bases = {}
                for pilups in p.pileups:
                    if pilups.query_position:
                        bp = pilups.alignment.query_sequence[pilups.query_position]
                        bases[bp] = bases.get(bp, 0) + 1
                # if we don't have any bases at this position, add an N
                if not bases:
                    bases['N'] = 1
                bps = sorted(bases, key=bases.get, reverse=True)
                # text = ""
                # for b in bps:
                #     text += " " + b + ": " + str(bases[b])
                # sys.stdout.write("{} : {} -> {}\n".format(p.reference_name, p.reference_pos, text))

                # make the consensus seq
                cseq.append(bps[0])
            print(">{} {} {} {}\n{}".format(primer, template, start, end, ''.join(cseq)))


def print_alignment(bam):
    """
    Print an alignment of all matching pileups

    :param bam: the bam object from pysam
    :type bam: pysam.AlignmentFile
    :return:
    :rtype:
    """

    for template in locations:
        for primer in locations[template]:
            start, end = locations[template][primer]
            alignment = {}
            for p in bam.pileup(reference=template, start=start, end=end, truncate=True):
                for pilups in p.pileups:
                    if pilups.alignment.query_name not in alignment:
                        alignment[pilups.alignment.query_name] = ['-' for idx in range(start, end)]
            for p in bam.pileup(reference=template, start=start, end=end, truncate=True):
                for pilups in p.pileups:
                    if pilups.query_position:
                        posn = pilups.query_position - start
                        alignment[pilups.alignment.query_name][posn] = \
                            pilups.alignment.query_sequence[pilups.query_position]

            # find the longest name
            longest_name = 0
            for n in alignment:
                if len(n) > longest_name:
                    longest_name = len(n)
            longest_name += 5

            for n in alignment:
                sys.stdout.write(n)
                sys.stdout.write(" " * (longest_name - len(n)))
                sys.stdout.write(''.join(alignment[n]) + "\n")

        print("\n\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract PCRd regions from BAM files')
    parser.add_argument('-b', help='bam file', required=True)
    parser.add_argument('-q', help='print query regions', action='store_true')
    parser.add_argument('-p', help='print pileup', action='store_true')
    parser.add_argument('-c', help='print consensus sequence', action='store_true')
    parser.add_argument('-a', help='print alignment', action='store_true')
    parser.add_argument('-v', help='verbose output')
    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.b, 'rb')

    if args.q:
        print_query_regions(bam)
    if args.p:
        print_pileup(bam)
    if args.c:
        print_consensus(bam)
    if args.a:
        print_alignment(bam)