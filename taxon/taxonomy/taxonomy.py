"""
Classes of taxonomy data
"""

import os
import sys


class NoNameFoundError(Exception):
    """No name was found for this entry"""
    def __init__(self, message):
        self.message = message


class TaxonNode:
    def __init__(self, t=None, p=None, r=None, e=None, d=None, i=None, gc=None, igc=False, mgc=None, imgc=False,
                 gh=False, hs=False, c=None, *others):
        self.parent = p
        self.taxid = t
        self.rank = r
        self.embl = e
        self.division = d
        self.inherited = i
        self.geneticCode = gc
        self.inheritedGC = igc
        self.mitochondrialGeneticCode = mgc
        self.inheritedMitochondrialGeneticCode = imgc
        self.GenBankHidden = gh
        self.hiddenSubtree = hs
        self.comments = c
        if len(others) > 0:
            print("WARNING: {} :: {}".format(p, others))


class TaxonName:
    def __init__(self, t=None, n=None, u=None, nc=None):

        self.taxid = t
        self.acronym = []
        self.anamorph = []
        self.authority = []
        self.blast_name = None
        self.common_name = []
        self.equivalent_name = []
        self.genbank_acronym = None
        self.genbank_anamorph = None
        self.genbank_common_name = None
        self.genbank_synonym = []
        self.in_part = []
        self.includes = []
        self.misnomer = []
        self.misspelling = []
        self.scientific_name = None
        self.synonym = []
        self.telemorph = []
        self.type_material = []
        self.unique = u

    def set_name(self, nametype, nameval):
        """
        Set the name for this tax id. ALlows multiple names for same ID
        :param nametype: the type of name
        :param nameval: the name
        """
        if "acronym" == nametype:
            self.acronym.append(nameval)
        elif "anamorph" == nametype:
            self.anamorph.append(nameval)
        elif "authority" == nametype:
            self.authority.append(nameval)
        elif "blast name" == nametype:
            self.blast_name = nameval
        elif "common name" == nametype:
            self.common_name.append(nameval)
        elif "equivalent name" == nametype:
            self.equivalent_name.append(nameval)
        elif "genbank acronym" == nametype:
            self.genbank_acronym = nameval
        elif "genbank anamorph" == nametype:
            self.genbank_anamorph = nameval
        elif "genbank common name" == nametype:
            self.genbank_common_name = nameval
        elif "genbank synonym" == nametype:
            self.genbank_synonym.append(nameval)
        elif "in-part" == nametype:
            self.in_part.append(nameval)
        elif "includes" == nametype:
            self.includes.append(nameval)
        elif "misnomer" == nametype:
            self.misnomer.append(nameval)
        elif "misspelling" == nametype:
            self.misspelling.append(nameval)
        elif "scientific name" == nametype:
            self.scientific_name = nameval
        elif "synonym" == nametype:
            self.synonym.append(nameval)
        elif "teleomorph" == nametype:
            self.telemorph.append(nameval)
        elif "type material" == nametype:
            self.type_material.append(nameval)
        else:
            sys.stderr.write("Do not recognise name type |{}|\n".format(nametype))

    def get_name(self):
        """
        Get the preferred name for this taxon
        :return: a string with the name
        """

        if self.blast_name:
            return self.blast_name
        if self.scientific_name:
            return self.scientific_name
        if self.common_name:
            return self.common_name
        if self.equivalent_name:
            return self.equivalent_name
        raise NoNameFoundError(f"No name was found for taxonomy ID {self.taxid}")

class TaxonDivision:
    def __init__(self, i=None, c=None, n=None, co=None):
        self.divid = i
        self.name = n
        self.code = c
        self.comments = co

