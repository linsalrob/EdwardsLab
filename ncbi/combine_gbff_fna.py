import os
import sys
import re
import roblib
import gzip

'''Combine .gbff and .fna files to get just the coding sequences. We need to get the data from RefSeq and they have 
split DNA sequences out of GenBank files so it is not clear that biopython etc will work.

This is just a quick parser and then we get the strings.'''


try:
    gbff = sys.argv[1]
    fnaf = sys.argv[2]
except:
    sys.stderr.write(sys.argv[0] + " <gbff file> <fna file>\n")
    sys.exit(-1)


locusre = re.compile('LOCUS\s+(\S+)')
locustagre = re.compile('\s+\/locus_tag=\"(.*)\"')
locationre = re.compile('\s+gene\s+(\d+)\.\.(\d+)$')
locationrerc = re.compile('\s+gene\s+complement\((\d+)\.\.(\d+)\)$')

locus = ""
locustag = ""
[start, end]=['0','0']
complement = False

locations={}

try:
    if gbff.endswith('.gz'):
        gbfin=gzip.open(gbff, 'rb')
    else:
        gbfin=open(gbff, 'r')
except:
    sys.exit("Unable to open file " + gbff)


for line in gbfin:
    line = line.rstrip()
    if line == "//":
        if start != '0' or end != '0':
            # print "\t".join([locus, locustag, start, end, str(complement)])
            locations[locus][locustag]=[start, end, complement]
        locus = ""
        locustag = ""
        [start, end]=['0','0']
        complement = False
        continue

    if line.startswith('LOCUS'):
        m = locusre.match(line)
        locus = m.group(1)
        locations[locus]={}
        continue

    if '/locus_tag' in line:
        m = locustagre.match(line)
        if m:
            locustag = m.group(1)
        else:
            sys.stderr.write("Couldn't parse |" + line + "|\n")

    if '..' in line and 'gene' in line:
        if start != '0' or end != '0':
            # print "\t".join([locus, locustag, start, end, str(complement)])
            locations[locus][locustag]=[start, end, complement]
        locustag = ""
        [start, end]=['0','0']
        complement = False

        m = locationre.match(line)
        if m:
            start = m.group(1)
            end = m.group(2)
        else:
            m = locationrerc.match(line)
            if m:
                complement = True
                start = m.group(1)
                end = m.group(2)
            else:
                sys.stderr.write("Can't parse an apparent location at : " + line + "\n")

fa = roblib.readFasta(fnaf)

ncre = re.compile('.*ref\|(\w+)')
for id in fa:
    m = ncre.match(id)
    if not m:
        sys.stderr.write("No apparent NC_ idenitifer in this sequence id: " + id + "\n")
        continue

    locus = m.group(1)
    for l in locations[locus]:
        [start, end, complement] = locations[locus][l]
        if complement:
            print ">" + l + " "  + locus + " " + end + "_" + start + " COMPLEMENT"
            print roblib.rc(fa[id][int(start) - 1:int(end)])
        else:
            print ">" + l + " "  + locus + " " + start + "_" + end
            print fa[id][int(start)-1:int(end)]



