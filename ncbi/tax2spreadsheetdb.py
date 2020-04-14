"""
Print the NCBI taxonomy as a spreadsheet
"""

from taxon import get_taxonomy_db, get_taxonomy

want = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']

def printtaxa(t, n, i, c):
    """
    Print out the taxonomy
    :param t: taxonomy object
    :param n: name
    :param i: identifier
    :param c: database connection
    :return:
    """

    names = {w: "" for w in want}

    names[t.rank] = n

    while t.parent != 1 and t.taxid != 1:
        t, n = get_taxonomy(t.parent, c)
        if t.rank in want:
            names[t.rank] = n
    print("\t".join([str(i)] + [names[w] for w in want]))



if __name__ == '__main__':
    c = get_taxonomy_db()
    for i in ids:
        t, n = get_taxonomy(i, c)
        if t.rank == "species":
            printtaxa(i)

