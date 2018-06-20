import sys
from taxon import get_taxonomy_db, get_taxonomy

## run this as a script
if __name__ == "__main__":
    ids=['1211997']
    if len(sys.argv) == 1:
        print("Supply your own id for the tree. Here is an example")
        print("""
              1211997 name: Cimicodes sp. BOLD:AAW5588    blast name: Cimicodes sp. BOLD:AAW5588   rank: species   parent: 704296
              704296 name: Cimicodes    blast name: Cimicodes   rank: genus   parent: 82596
              82596 name: Ennominae    blast name: Ennominae   rank: subfamily   parent: 82593
              82593 name: Geometridae    blast name: Geometridae   rank: family   parent: 82592
              82592 name: Geometroidea    blast name: Geometroidea   rank: superfamily   parent: 104431
              104431 name: Obtectomera    blast name: Obtectomera   rank: no rank   parent: 37567
              37567 name: Ditrysia    blast name: Ditrysia   rank: no rank   parent: 41197
              41197 name: Heteroneura    blast name: Heteroneura   rank: parvorder   parent: 41196
              41196 name: Neolepidoptera    blast name: Neolepidoptera   rank: infraorder   parent: 41191
              41191 name: Glossata    blast name: Glossata   rank: suborder   parent: 7088
              7088 name: Lepidoptera    blast name: moths   rank: order   parent: 85604
              85604 name: Amphiesmenoptera    blast name: Amphiesmenoptera   rank: no rank   parent: 33392
              33392 name: Holometabola    blast name: Holometabola   rank: cohort   parent: 33340
              33340 name: Neoptera    blast name: Neoptera   rank: infraclass   parent: 7496
              7496 name: Pterygota    blast name: Pterygota   rank: subclass   parent: 85512
              85512 name: Dicondylia    blast name: Dicondylia   rank: no rank   parent: 50557
              50557 name: Insecta    blast name: Insecta   rank: class   parent: 6960
              6960 name: Hexapoda    blast name: insects   rank: superclass   parent: 197562
              197562 name: Pancrustacea    blast name: Pancrustacea   rank: no rank   parent: 197563
              197563 name: Mandibulata    blast name: Mandibulata   rank: no rank   parent: 6656
              6656 name: Arthropoda    blast name: arthropods   rank: phylum   parent: 88770
              88770 name: Panarthropoda    blast name: Panarthropoda   rank: no rank   parent: 1206794
              1206794 name: Ecdysozoa    blast name: Ecdysozoa   rank: no rank   parent: 33317
              33317 name: Protostomia    blast name: Protostomia   rank: no rank   parent: 33213
              33213 name: Bilateria    blast name: Bilateria   rank: no rank   parent: 6072
              6072 name: Eumetazoa    blast name: Eumetazoa   rank: no rank   parent: 33208
              33208 name: Metazoa    blast name: animals   rank: kingdom   parent: 33154
              33154 name: Opisthokonta    blast name: Opisthokonta   rank: no rank   parent: 2759
              2759 name: Eukaryota    blast name: eukaryotes   rank: superkingdom   parent: 131567
              """)
        sys.exit(0)
    else:
        ids=sys.argv[1:]


    c = get_taxonomy_db()
    for i in ids:
        t, n = get_taxonomy(i, c)
        while t.parent != 1 and t.taxid != 1:
            print("{}\tname: {}\trank: {}\tparent: {}".format(t.taxid, n.scientific_name, t.rank, t.parent))
            t, n = get_taxonomy(t.parent, c)


