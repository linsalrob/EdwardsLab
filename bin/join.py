import os
import sys

__author__ = 'Rob Edwards'

data={}
allhs = set()
for f in sys.argv[1:]:
    sys.stderr.write("Parsing {}\n".format(f))
    with open(f, 'r') as fn:
        header = None
        for l in fn:
            if not header:
                header=l.strip().split("\t")
                allhs.update(set(header))
                continue
            p = l.strip().split("\t")
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
        sys.stdout.write("\t{}".format(data[seqlib].get(contig, 0)))
    sys.stdout.write("\n")