import sys

sl = [] # sorted list
with open("RecA_uniprot_pairwise_ids.tsv", 'r') as f:
#with open("temp", 'r') as f:
    for l in f:
        p=l.strip().split("\t")
        val = float(p[2])
        added = False
        for i in range(len(sl)):
            if sl[i] > val:
                sl[i:i]=[val]
                added=True
                break
        if not added:
            sl.append(val)

indices = [-1 for i in range(100)]
for i in range(len(sl)):
    ni = int(sl[i])+1
    sys.stderr.write("Adding {} for {}\n".format(ni, sl[i]))
    indices[ni]=i

print(sl)

lastindex=-1
for i in range(50, 70):
    if indices[i] == -1:
        indices[i] = lastindex
    else:
        lastindex = indices[i]
    print("{}\t{}".format(i, indices[i]+1))

