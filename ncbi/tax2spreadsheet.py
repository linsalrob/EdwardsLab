import taxon

taxa=taxon.readNodes()
names,blastname = taxon.readNames()
divs = taxon.readDivisions()

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

    print i  + "\t" + bn,
    for l in want:
        if l in level:
            print "\t" + level[l],
        else:
            print "\t-",
    print



# levels: species genus family order class phylum  kingdom

print "id\tname",
for l in want:
    print "\t" + l,
print

for i in taxa:
    if taxa[i].rank == "species":
        printtaxa(i)

