import argparse
import sys

import pysam

__author__ = 'Rob Edwards'

"""
Read a BAM file and calculate the coverage for each read
"""






if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='calculate the coverage of each read in a bam file')
    parser.add_argument('-b', help='bam file to parse')
    parser.add_argument('-s', help='sam file to parse')
    args = parser.parse_args()

    samfile = None
    if args.b:
        samfile = pysam.AlignmentFile(args.b, "rb")
    elif args.s:
        samfile = pysam.AlignmentFile(args.s, "r")
    else:
        sys.exit("Either -s or -b must be specified")

    alllibs = set()
    count = {}
    for read in samfile.fetch():
        #     print("{} -> {} : {}".format(read.query_name, read.reference_name, read.query_alignment_length))
        if '#' not in read.query_name:
            sys.stderr.write("# not found in {}\n".format(read.query_name))
            continue

        library =  read.query_name.split('#')[1].split('/')[0]

        if read.reference_name not in count:
            count[read.reference_name] = {}
        if library not in count[read.reference_name]:
            count[read.reference_name] = {}

        count[read.reference_name][library] = count[read.reference_name].get(library, 0)+1

        alllibs.add(library)

    als = list(alllibs)
    als.sort()

    print("\t" + "\t".join(als))
    for r in count.keys():
        sys.stdout.write(r)
        for l in als:
            sys.stdout.write("\t{}".format(count[r].get(l, 0)))
        sys.stdout.write("\n")

