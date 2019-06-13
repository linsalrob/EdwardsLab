from taxon import get_taxonomy_db, get_taxonomy, all_ids

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


c = get_taxonomy_db()
for i in all_ids(c):
    print (f"{i}")
    t, n = get_taxonomy(i, c)
    if t.rank == "phylum":
        while t.parent != 1 and t.taxid != 1:
            t, n = get_taxonomy(t.parent, c)
            print(f"rank: {t.rank} :: name: {t.common_name}")


