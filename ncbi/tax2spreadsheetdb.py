"""
Print the NCBI taxonomy as a spreadsheet
"""

from taxon import get_taxonomy_db, get_taxonomy, all_species_ids

want = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

def printtaxa(i, c):
    """
    Print out the taxonomy
    :param i: identifier
    :param c: database connection
    :return:
    """

    names = {w: "" for w in want}
    t, n = get_taxonomy(i, c)
    if t.rank in want:
        names[t.rank] = n
    while t.parent != 1 and t.taxid != 1:
        t, n = get_taxonomy(t.parent, c)
        if t.rank in want:
            names[t.rank] = n.get_name()
    print("\t".join([str(i)] + [names[w] for w in want]))



if __name__ == '__main__':
    c = get_taxonomy_db()
    for i in all_species_ids(c):
        printtaxa(i[0], c)
