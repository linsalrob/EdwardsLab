import taxon

taxa=taxon.read_nodes()
names,blastname = taxon.read_names()
divs = taxon.read_divisions()

want = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']

def printtaxa(i):
    bn=names[i].name
    if i in blastname:
        bn=blastname[i].name
    
    level={}
    
    node = i
    while taxa[node].parent != '1' and node != '1':
        if taxa[node].rank in want:
            level[taxa[node].rank]=names[node].name
        node=taxa[node].parent

    print("{}\t{}".format(i, bn), end="")
    for l in want:
        if l in level:
            print("\t{}".format(level[l]), end="")
        else:
            print("\t-", end="")
    print("")



# levels: species genus family order class phylum  kingdom

print ("id\tname", end="")
for l in want:
    print("\t{}".format(l), end="")
print("")

for i in taxa:
    if taxa[i].rank == "species":
        printtaxa(i)

