import argparse
import os
import sys
import pysam

__author__ = 'Rob Edwards'

"""
This is just a version of extract-pcr-regions adapted to print out information from Ondrej Cinek.
Mainly it is the consensus header. You can probably delete this file!

"""

locations = {
    'JQ995537': {
        'A': (25634, 26964),
        'B': (33709, 35062),
        'C': (43820, 45057)
    }
}


def print_consensus(bam, sname):
    """
    Print the consensus sequence about the pileup

    :param bam: the bam object from pysam
    :type bam: bam
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
            header = 'Cinek_{}_primer{}.20151120 [name=Cinek lab sample {} primer {} 20151120] '.format(sname, primer, sname, primer)
            header += '[date=20151120] [latitude=61.505] [longitude=23.815] [note=coordinates are hospital location]'
            if cseq:
                print(">{}\n{}".format(header, ''.join(cseq)))

# >Brouns_A.20151106 [name=Brouns lab primer A 20151106] [date=20151106] [latitude=52.180646] [longitude=5.9478529] [altitude=36m]



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract PCRd regions from BAM files')
    parser.add_argument('-b', help='bam file', required=True)
    parser.add_argument('-v', help='verbose output')
    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.b, 'rb')

    sname = args.b.split('.')[0]

    # print_query_regions(bam)
    # print_pileup(bam)
    print_consensus(bam, sname)