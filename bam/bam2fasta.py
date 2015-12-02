import argparse

import pysam

__author__ = 'Rob Edwards'


def qual2fastq(quals):
    """
    Convert a list of quality scores to a single fastq line

    :param quals: A list of quality scores
    :type quals: list
    :return: A fastq quality string
    :rtype: str
    """
    quality = [chr(q) for q in quals]
    return "".join(quality)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert bam to fastq')
    parser.add_argument('-b', help='bam file', required=True)
    parser.add_argument('-v', help='verbose output')
    args = parser.parse_args()

    bamfile = pysam.AlignmentFile(args.b, "rb")
    for read in bamfile.fetch(until_eof=True):
        print(">{}\n{}".format(read.query_name, read.query_sequence))
