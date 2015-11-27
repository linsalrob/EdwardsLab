import os
import sys

__author__ = 'Rob Edwards'

data={}
allhs = set()
for f in sys.argv[1:]:
    sys.stderr.write("Parsing {}\n".format(f))
    with open(f, 'r') as fn:
        header = []
        for l in fn:
            if header == []:
                header=l.strip().split("\t")
                allhs.update(set(header))
                header.insert(0, "Contig")
                sys.stderr.write("Set header: {}\n".format(header))
                continue
            p = l.strip().split("\t")
            if len(p) != len(header):
                sys.exit("P has length: {} and header has length {}".format(len(p), len(header)))
            if p[0] not in data:
                data[p[0]] = {}
            for i in range(1, len(p)):
                data[p[0]][header[i]] = data[p[0]].get(header[i], 0) + int(p[i])

allhsl = list(allhs)
allhsl.sort()
print("\t" + "\t".join(allhsl))
for contig in data:
    sys.stdout.write(contig)
    for seqlib in allhsl:
        sys.stdout.write("\t{}".format(data[contig].get(seqlib, 0)))
    sys.stdout.write("\n")
