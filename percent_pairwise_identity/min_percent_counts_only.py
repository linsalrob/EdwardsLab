"""

This code is designed to work with RecA_percents_sort.txt and integrase_percents_sort.txt

These are text files that have one entry per line which is the percent ID and are forward sorted (so the highest number
is last). They are derived from RecA_uniprot_pairwise_ids.tsv and
integrase.similarity.tsv by using  cut -f 3 -d$'\t' integrase.similarity.tsv | sort -n

The percent IDs are already sorted, and there is one entry per line.

"""



import sys

sl = [] # sorted list
# 220000

counter = 0; oco = 0;
#with open("RecA_percents_sort.txt", 'r') as f:
with open("integrase_percent_sort.txt", 'r') as f:
    sys.stderr.write("File open. Lets go\n")
    for l in f:
        counter += 1
        if (counter % 220000) == 0:
            oco += 1
            sys.stderr.write("Read {} : {}\n".format(oco, counter))
        val = float(l.strip())
        sl.append(val)

indices = [-1 for i in range(110)]
for i in range(len(sl)):
    ni = int(sl[i])+1
    indices[ni]=i

lastindex=-1
for i in range(1, 102):
    if indices[i] == -1:
        indices[i] = lastindex
    else:
        lastindex = indices[i]
    print("{}\t{}".format(i-1, indices[i]+1))

