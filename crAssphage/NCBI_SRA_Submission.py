"""
Submit my data to the SRA.
"""

import os
import sys
import argparse
from roblib import bcolors
__author__ = 'Rob Edwards'



our_data = {
    'library_strategy'    : 'AMPLICON',
    'library_source'      : 'METAGENOMIC',
    'library_selection'   : 'PCR',
    'library_layout'      : 'SINGLE',
    'platform'            : 'capillary',
    'instrument_model'    : 'AB 5500 Genetic Analyzer',
}


def read_sra_data(ifile='sra_ids_experiments.txt', verbose=False):
    """
    Read the SRA experiment ID that we generated with an SQL command
    """

    sra_data = {}
    # library_strategy,library_source,library_selection,library_layout,platform,instrument_model,platform_parameters
    with open(ifile, 'r') as f:
        for l in f:
            p = l.rstrip("\n").split("\t")
            sra_data[p[0]] = {}
            for i,j in enumerate(['library_strategy', 'library_source', 'library_selection', 'library_layout', 'platform', 'instrument_model', 'platform_parameters']):
                sra_data[p[0]][j]=p[i+1] # the first element is the id :)
    return sra_data


def parse_file(ifile, ofile, sra_data, verbose=False):
    """
    Parse the file and create a new metadata file
    """
    cols = {'bioproject_accession' : 0, 'biosample_accession' : 1, 'library_ID' : 2, 
            'title' : 3, 'library_strategy' : 4, 'library_source' : 5, 'library_selection' : 6,
            'library_layout' : 7, 'platform' : 8, 'instrument_model' : 9, 'design_description' : 10,
            'filetype' : 11, 'filename' : 12, 'filename2' : 13, 
            'filename3' : 14, 'filename4' : 15, 'assembly' : 16}
    
    with open(ifile, 'r') as f:
        with open(ofile, 'w') as out:
            out.write("bioproject_accession\tbiosample_accession\tlibrary_ID\ttitle\tlibrary_strategy\tlibrary_source\tlibrary_selection\tlibrary_layout\tplatform\tinstrument_model\tdesign_description\tfiletype\tfilename\tfilename2\tfilename3\tfilename4\tassembly\n")
            for l in f:
                if l.startswith("BioProject_Accession"):
                    continue
                p = l.rstrip("\n").split("\t")
                # important columns:
                    # bioproject_accession: 0
                    # biosample_accession: 1
                    # library_ID: 2
                    # title: 188
                    # src: 200
                    # sra_id: 199

                row = [p[0], p[1], p[2], p[188], "", "", "", "", "", "", "", "fastq", "crAssphage.fasta", "", "", "", ""]

                if 'SRA' == p[200]:
                    if not p[199]:
                        sys.stderr.write(f"{bcolors.FAIL}FATAL: for {p[2]} src is SRA but there is no SRA ID{bcolors.ENDC}\n")
                        sys.exit(-1)
                    if p[199] not in sra_data:
                        sys.stderr.write(f"{bcolors.FAIL}FATAL: for {p[2]} SRA ID {p[199]} is not in the sra data{bcolors.ENDC}\n")
                        sys.exit(-1)
                    for k in ['library_strategy', 'library_source', 'library_selection', 'library_layout', 'platform', 'instrument_model']:
                        if k not in sra_data[p[199]]:
                            sys.stderr.write(f"{bcolors.FAIL}FATAL: {k} not in SRA data for {p[199]}\n")
                            continue
                        row[cols[k]] = sra_data[p[199]][k]
                else:
                    for k in our_data:
                        row[cols[k]] = our_data[k]

                out.write("{}\n".format("\t".join(row)))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='all_data biosamples file', required=True)
    parser.add_argument('-o', help='file to write', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    
    sra_data = read_sra_data('sra_ids_experiments.txt', args.v)
    parse_file(args.f, args.o, sra_data, args.v)


