"""
Convert the id.map file and seqs.dnadist file to a single matrix that we can use in an anova to test
the relationship between sequences and locations

Our current keys in the metadata are:
    address
    altitude
    country
    date
    latitude
    longitude
    name
    note
    source

We want to use name, country, date, latitude, longitude as our keys in the output


"""

import os
import sys
import argparse
import re

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert seqs.dnadist and seqs.')
    parser.add_argument('-d', help='seqs.dnadist file', required=True)
    parser.add_argument('-i', help='id.map file', required=True)
    parser.add_argument('-o', help='output file', required=True)
    args = parser.parse_args()

    metadata={}
    allmd=set()
    allmd.add('seqid')

    # read the id.map file
    with open(args.i, 'r') as f:
        for l in f:
            p=l.strip().split("\t")
            metadata[p[0]] = {}

            m=re.match('[^\[]*', p[1])
            seqid = m.group().strip()
            if not seqid:
                sys.stderr.write("Can't parse a seqid from {}\n".format(l))
                sys.exit()

            metadata[p[0]]['seqid'] = seqid

            m=re.findall('\[.*?\]', p[1])
            for s in m:
                s = s.replace('[', '')
                s = s.replace(']', '')
                d = s.split('=')
                metadata[p[0]][d[0]]=d[1]
                allmd.add(d[0])


    dnadists = []
    thisdist = ""
    first = True
    with open(args.d, 'r') as f:
        for l in f:
            if first:
                # ignore the first line
                first = False
                continue
            if not l.startswith(' '):
                if thisdist:
                    dnadists.append(thisdist)
                thisdist = l.strip()
            else:
                thisdist += ' ' + l.strip()
    dnadists.append(thisdist)

    with open(args.o, 'w') as out:
        for l in dnadists:
                p=l.strip().split()
                if '_R_' in p[0]:
                    p[0] = p[0].replace('_R_', '')
                if p[0] not in metadata:
                    sys.stderr.write("ERROR: No metadata for {} in {}\n".format(p[0], l))
                    sys.exit()

                for tag in ['name', 'date', 'latitude', 'longitude', 'country']:
                    if tag in metadata[p[0]]:
                        out.write("{}\t".format(metadata[p[0]][tag]))
                    else:
                        out.write("\t")
                out.write("\t".join(map(str, p[1:])))
                out.write("\n")





