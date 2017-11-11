import os
import sys
import taxon
import argparse


taxa = taxon.read_nodes()
names, blastname = taxon.read_names()
divs = taxon.read_divisions()

indices = {names[i].name.lower():i for i in names}


def plant_or_animal(namestr):
    """
    Determine whether something is a plant or an animal or an insect
    :param namestr:
    :return:
    """

    namestr = namestr.lower()

    if namestr not in indices:
        return None

    i = indices[namestr]

    while taxa[i].parent != '1' and i != '1':
        if names[i].name in ['Metazoa', 'Viridiplantae', 'Arthropoda']:
            return names[i].name
        i = taxa[i].parent



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Try and associate plants and animals with names")
    parser.add_argument('-f', help='file of names to search for')
    args = parser.parse_args()

    with open(args.f, 'r') as f:
        for l in f:
            n = plant_or_animal(l.strip())
            if not n:
                n = plant_or_animal(l.split(' ')[0].replace(',', ''))
            print("{}\t{}".format(l.strip(), n))











