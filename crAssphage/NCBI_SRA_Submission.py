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
    'instrument_model'    : 'AB 3730 Genetic Analyzer'
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
            if 'PAIRED' in p[4]:
                p[4] = 'PAIRED'
            if 'SINGLE' in p[4]:
                p[4] = 'SINGLE'
            if 'unspecified' in p[6]:
                p[6] = 'Illumina HiSeq 1000'
            for i,j in enumerate(['library_strategy', 'library_source', 'library_selection', 'library_layout', 'platform', 'instrument_model', 'platform_parameters']):
                sra_data[p[0]][j]=p[i+1] # the first element is the id :)
    return sra_data


def parse_file(ifile, ofile, sra_data, out_dir, verbose=False):
    """
    Parse the file and create a new metadata file.

    Writes everything to the directory SRA_Submission, including a set of fasta files, 
    one per biosample.

    """
    cols = {'bioproject_accession' : 0, 'biosample_accession' : 1, 'library_ID' : 2, 
            'title' : 3, 'library_strategy' : 4, 'library_source' : 5, 'library_selection' : 6,
            'library_layout' : 7, 'platform' : 8, 'instrument_model' : 9, 'design_description' : 10,
            'filetype' : 11, 'filename' : 12, 'filename2' : 13, 
            'filename3' : 14, 'filename4' : 15, 'assembly' : 16}
    

    if os.path.exists(out_dir):
        sys.stderr.write(f"{bcolors.FAIL}ERROR: {out_dir} exists. Not overwriting\n")
        sys.exit(-1)

    os.mkdir(out_dir)

    volume = {}

    linecount = 1
    filecount = 0
    if ".tsv" in ofile:
        ofile = ofile.replace(".tsv", "")
    
    out = open(f"{out_dir}/{ofile}.{filecount}.tsv", 'w')
    out.write("bioproject_accession\tbiosample_accession\tlibrary_ID\ttitle\tlibrary_strategy\tlibrary_source\tlibrary_selection\tlibrary_layout\tplatform\tinstrument_model\tdesign_description\tfiletype\tfilename\tfilename2\tfilename3\tfilename4\tassembly\n")
    os.mkdir(os.path.join(out_dir, str(filecount)))
    increment_filecount = False

    with open(ifile, 'r') as f:
        for l in f:
            if l.startswith("BioProject_Accession"):
                continue

            # this is to ensure that we write the fasta sequences to the correct subdir
            # for sequences that we are not processing further
            if increment_filecount:
                filecount+=1
                increment_filecount = False
                out.close()
                out = open(f"{out_dir}/{ofile}.{filecount}.tsv", 'w')
                out.write("bioproject_accession\tbiosample_accession\tlibrary_ID\ttitle\tlibrary_strategy\tlibrary_source\tlibrary_selection\tlibrary_layout\tplatform\tinstrument_model\tdesign_description\tfiletype\tfilename\tfilename2\tfilename3\tfilename4\tassembly\n")
                os.mkdir(os.path.join(out_dir, str(filecount)))


            p = l.rstrip("\n").split("\t")
            # important columns:
                # bioproject_accession: 0
                # biosample_accession: 1
                # sequence_ID: 2
                # sample name: 3
                # title: 188
                # src: 200
                # sra_id: 199
                # sequence: 234

            # write the sequence out
            # note that we do this before we process the line because we need 
            # all sequences, but we don't process all of them if we have already
            # seen the biosample ID before
            subdir = str(filecount)
            if p[1] in volume:
                subdir = str(volume[p[1]])
            fao = open(os.path.join(out_dir, subdir, f"{p[1]}.fasta"), 'a')
            fao.write(f">{p[2]}\n{p[234]}\n")
            fao.close()


            if p[1] in volume:
                continue
            volume[p[1]] = filecount

            linecount+=1
            if linecount > 999:
                linecount = 1
                increment_filecount = True
            
            row = [p[0], p[1], p[3], p[188], "", "", "", "", "", "", "", "fastq", f"{p[1]}.fasta", "", "", "", ""]

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
                row[cols['design_description']] = f"Extracted from SRA run {p[199]} using Gretel version 0.0.8 (https://github.com/SamStudio8/gretel)" 
            else:
                for k in our_data:
                    row[cols[k]] = our_data[k]
                row[cols['design_description']] = f"PCR amplified from a raw environmental sample and seqeunced using Sanger sequencing. "

            out.write("{}\n".format("\t".join(row)))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='all_data biosamples file', required=True)
    parser.add_argument('-o', help='file to write', required=True)
    parser.add_argument('-d', help='directory to write to', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    
    sra_data = read_sra_data('sra_ids_experiments.txt', args.v)
    parse_file(args.f, args.o, sra_data, args.d, args.v)


