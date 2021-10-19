"""
Convert a two dimensional matrix to an h5 format file for use with turbocor

This takes a tab, or comma, separated file and converts it into a 2D array
for use with turbocor.

Optionally, we write an index file with the row indices and can skip the
column headers.
"""

import os
import sys
import argparse
import numpy as np
import h5py

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse a tab or csv separated file and output it as an HDF5 format 2D array')
    parser.add_argument('-f', '--file', help='input tsv/csv file to convert', required=True)
    parser.add_argument('-o', '--output', help='output hdf5 format file to write', required=True)
    parser.add_argument('-r', '--header', help='Skip the first row because it is a header', action='store_true')
    parser.add_argument('-d', '--dataset', help='dataset name that will be used in the output (default = data)', default='data')
    parser.add_argument('-s', '--separator', help='use an alternate input record separator(default: tab)', default="\t")
    parser.add_argument('-n', '--npdtype', help='default data type for the input data (default=int). See numpy dtype for codes', default='i')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-c', '--col', help='Include the first column (default: use it as the ID)', action='store_true')
    group.add_argument('-i', '--indexfile', help='index file of first column IDs and position to write (optional but makes parsing output easier!)')

    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    dtypes = {
        'i': 'integer',
        'b': 'boolean',
        'u': 'unsigned integer',
        'f': 'float',
        'c': 'complex float',
        'm': 'timedelta',
        'M': 'datetime',
        'O': 'object',
        'S': 'string',
        'U': 'unicode string',
        'V': 'void'
    }

    if args.npdtype not in dtypes:
        sys.stderr.write(f"FATAL: npdtype must be one of {dtypes}")
        sys.exit(1)

    data = []
    header = False
    with open(args.file, 'r') as f:
        if args.indexfile:
            indexout = open(args.indexfile, 'w')
        counter = -1
        for l in f:
            if args.header and not header:
                header = True
                continue
            counter += 1
            p = l.strip().split(args.separator)
            if not args.col:
                colname = p.pop(0)
                if args.indexfile:
                    indexout.write(f"{counter}\t{colname}\n")
            data.append(np.array(p, dtype=args.npdtype))
        if args.indexfile:
            indexout.close()

    with h5py.File(args.output, "w") as f:
        grp = f.create_dataset(args.dataset, data=np.array(data))

