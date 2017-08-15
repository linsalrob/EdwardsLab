import sys
import taxon

## run this as a script
if __name__ == "__main__":
    ids=['1211997']
    if len(sys.argv) == 1:
        print("Supply your own id for the tree. Here is an example")
    else:
        ids=sys.argv[1:]

    taxa = taxon.read_nodes()
    names,blastname = taxon.read_names()
    divs = taxon.read_divisions()
    for i in ids:
        c=0
        if i not in taxa:
            sys.stderr.write("Error: ID " + str(i) + " was not found in the taxonomy filei\n")
            sys.exit(0)
        while taxa[i].parent != '1' and i != '1':
            #print " " * c, taxa[i].taxid, "name:", names[i].name, "rank:", taxa[i].rank, "div:", taxa[i].division, "code:", divs[taxa[i].division].code
            bn=names[i].name
            if i in blastname:
                bn=blastname[i].name
            print("{} name: {}    blast name: {}   rank: {}   parent: {}".format(taxa[i].taxid, names[i].name, bn, taxa[i].rank, taxa[i].parent))
            i=taxa[i].parent
            c+=1


