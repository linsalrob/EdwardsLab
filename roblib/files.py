"""
Helper functions for files
"""
import binascii


def is_gzip(filename: str) -> bool:
    """
    Is this a gzip file?
    """

    """
    This is an elegant solution to test whether a file is gzipped by reading the first two characters.
    I also use a version of this in fastq_pair if you want a C version :)
    See https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed for inspiration
    :param f: the file to test
    :return: True if the file is gzip compressed else false
    """
    with open(filename, 'rb') as i:
        return binascii.hexlify(i.read(2)) == b'1f8b'
