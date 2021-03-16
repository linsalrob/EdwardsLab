import sys
from taxon import get_taxonomy_db, get_taxonomy


want = ['taxonomy', 'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']




## run this as a script
if __name__ == "__main__":
    ids=['1211997']
    if len(sys.argv) == 1:
        print("Supply your own id for the tree. Here is an example")
        print("\t".join(want))
        print("1211997\tEukaryota\tMetazoa\tArthropoda\tInsecta\tLepidoptera\tGeometridae\tCimicodes\tCimicodes")
        sys.exit(0)
    else:
        ids=sys.argv[1:]

    ids = sorted(ids)
    c = get_taxonomy_db()
    for i in ids:
        results = [i, '-', '-', '-', '-', '-', '-', '-', '-', '-']
        t, n = get_taxonomy(i, c)
        # while t.parent != 1 and t.taxid != 1:
        while t.taxid != 1:
            if t.rank in want:
                results[want.index(t.rank)] = n.scientific_name
            t, n = get_taxonomy(t.parent, c)
        print("\t".join(map(str, results)))


