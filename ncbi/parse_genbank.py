"""
Parse a genbank entry that we download from NCBI or a file that is on our hard drive.

This example uses BioPython
"""

import os
import sys
from Bio import Entrez, SeqIO
import re
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Parse a GenBank record from NCBI or a local file. Either -g or -f are required")
    parser.add_argument('-g', help='GenBank ID to download. This should just be the GI number')
    parser.add_argument('-f', help='Local file')
    parser.add_argument('-o', help='Output file for record downloaded from GenBank. Optional. If you want to save the entry for later')
    args = parser.parse_args()

    # define seq here, but we don't set it until later
    seq = ""

    if args.g:
        # retrieve a GI number from GenBank
        Entrez.email = 'raedwards@gmail.com' # set this so NCBI knows who to complain to

        handle = Entrez.efetch(db="nucleotide", id=args.g, rettype="gbwithparts", retmode="text")
        # we either write it to a file or we just read the stream
        if args.o:
            with open(args.o, 'w') as out:
                out.write(handle.read())
                handle.close()
            seq = SeqIO.read(args.o, "genbank")
        else:
            seq = SeqIO.read(handle, "genbank")
            handle.close()
    elif args.f:
        seq = SeqIO.read(args.f, 'genbank')
    else:
        sys.exit('Please provide either a filename (-f) or a GI number (-g) to parse')

    # seq is now or Seq object. We just need to parse what we need.



    ## this block will print locus tag, gene symbol, and dbxref for all features
    ## Remove the """ at the start and end of the block to use it.
    """
    for feature in seq.features:
        lt = 'None'
        if 'locus_tag' in feature.qualifiers:
            # there should only be one locus tag so we just use the first element in the list
            lt = feature.qualifiers['locus_tag'][0]

        gs = "None"
        if 'gene_synonym' in feature.qualifiers:
            gs = '; '.join(feature.qualifiers['gene_synonym'])

        db = "None"
        if 'db_xref' in feature.qualifiers:
            db = '; '.join(feature.qualifiers['db_xref'])

        print("\t".join([lt, gs, db]))
    """

    ## This block will just print the locus_tag, ECK, and GeneID for all features

    locus_tag = {}

    for feature in seq.features:
        lt = 'None'
        if 'locus_tag' in feature.qualifiers:
            lt = feature.qualifiers['locus_tag'][0]
        else:
            # we don't want to print all the information
            continue

        gene = "None"
        if 'gene' in feature.qualifiers:
            gene = feature.qualifiers['gene'][0]

        gs = "None"
        if 'gene_synonym' in feature.qualifiers:
            for gene_symbols in feature.qualifiers['gene_synonym']:
                # gene_symbols is a list of gene symbols joined with '; '
                for g in gene_symbols.split('; '):
                    if g.startswith("ECK"):
                        gs = g

        db = "None"
        if 'db_xref' in feature.qualifiers:
            for d in feature.qualifiers['db_xref']:
                if d.startswith('GeneID'):
                    db = d

        if lt in locus_tag:
            if gene != "None":
                locus_tag[lt][0] = gene
            if gs != "None":
                locus_tag[lt][1] = gs
            if db != "None":
                locus_tag[lt][2] = db
        else:
            locus_tag[lt] = [gene, gs, db]


    for l in locus_tag:
        print(l + "\t" + "\t".join(locus_tag[l]))
#        print("\t".join([lt, gene, gs, db]))

