"""
Download GenBank files from the GenBank assembly archive, and then parse
out the prophage and non-prophage regions.

You will need:
    - The PhiSpy locations file from https://open.flinders.edu.au/ndownloader/files/39664894  
    - The GenBank Assembly Summary file from https://open.flinders.edu.au/ndownloader/files/39664879

To get those, use:

    curl -Lo phage_locations_202206201.tsv.gz https://open.flinders.edu.au/ndownloader/files/39664894
    curl -Lo assembly_summary_20220601.txt.gz https://open.flinders.edu.au/ndownloader/files/39664879

Note:
    These are the versions of the files that we used in PhiSpy. You can
    also download a newer version of the GenBank Assembly Summary File
    from https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/


To run this script, use:

python download_phage_slices.py -f phage_locations_202206201.tsv.gz -a assembly_summary_20220601.txt.gz -p phage_regions -b bacteria_regions

You can also limit the number of phage contigs with the -m parameter.
The assembly summary file does not list the total number of contigs in the
assembly, and I have not included that information in these two files,
so we are limited to just the number of contigs on which we identified
potential phage regions. This is a proxy for the total number of
contigs in the genome. (Ask Rob for that info, if you really want it).

The -v flag will give you some more progress on the downloads, as it
lists each genome being processed.

Rob Edwards
07/07/23
Sydney Airport
"""

import errno
import os
import sys
import argparse
from io import StringIO
import requests
import gzip
import zlib
from Bio import SeqIO
__author__ = 'Rob Edwards'





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', '--phispy', help='phispy File summary of genomes, contigs and phage locations', required=True)
    parser.add_argument('-a', '--assembly', help='genbank assembly summary file', required=True)
    parser.add_argument('-p', '--phages', help='directory to write the phages to', required=True)
    parser.add_argument('-b', '--bacteria', help='directory to write the bacteria to', required=True)
    parser.add_argument('-m', '--maxcontigs', help='maximum number of contigs in this genome [Default: %(default)d]', default=5, type=int)
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()


    # read the phispy file
    if not os.path.exists(args.phispy):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.phispy)

    if args.verbose:
        print("Reading Phages", file=sys.stderr)
    phages = {}
    contig_count = {}
    with gzip.open(args.phispy, "rt") as f:
        tphages = {}
        for l in f:
            (genome, contig, start, stop, length, numgenes, reason) = l.strip().split("\t")
            if genome not in contig_count:
                contig_count[genome] = set()
            
            contig_count[genome].add(contig)

            if reason == "Kept":
                if genome not in tphages:
                    tphages[genome] = {}
                if contig not in tphages[genome]:
                    tphages[genome][contig] = []
                tphages[genome][contig].append([int(start), int(stop)])
        
        phages = {g:tphages[g] for g in tphages if len(contig_count[g]) <= args.maxcontigs}

    print(f"There are {len(phages.keys())} genomes with fewer than {args.maxcontigs} phage containing contigs", file=sys.stderr)

    #load the urls
    if args.verbose:
        print("Reading URLs", file=sys.stderr)
    urls = {}
    with gzip.open(args.assembly, "rt") as f:
        for l in f:
            p = l.strip().split("\t")
            if p[0] in phages:
                urls[p[0]] = p[19]


    os.makedirs(args.bacteria, exist_ok=True)
    os.makedirs(args.phages, exist_ok=True)

    for i,g in enumerate(phages.keys()):
        if args.verbose:
            print(f"Downloading phage {i}:  {g}", file=sys.stderr)
        asm = urls[g].split("/")[-1]
        url = f"{urls[g]}/{asm}_genomic.gbff.gz"
        r = requests.get(url)
        if r.status_code != 200:
            print(f"ERROR: Could not download {url}. Request response code was {r.status_code}", file=sys.stderr)
            continue
        bacteria_filehandle = open(os.path.join(args.bacteria, f"{g}.gbk"), "w")
        phage_filehandle = open(os.path.join(args.phages, f"{g}.gbk"), "w")
        with StringIO(zlib.decompress(r.content, 16+zlib.MAX_WBITS).decode('utf-8')) as gio:
            for record in SeqIO.parse(gio, "genbank"):
                if record.id in phages[g]:
                    hostcounter = 0
                    ppnum = 0
                    bactstart = 0
                    bactend = -1
                    phagestart = -1
                    phageend = -1
                    for loc in phages[g][record.id]:
                        bactend = loc[0] - 1
                        phagestart = loc[0]
                        phageend = loc[1]

                        host_gbk = record[bactstart:bactend]
                        hostcounter += 1
                        host_gbk.name += f"_region_{hostcounter}"
                        host_gbk.id += f"_region_{hostcounter}"
                        host_gbk.description += f" region {hostcounter} on {contig} from {bactstart} to {bactend}"                                                               
                        # set the molecule type; required for Biopython v1.78 and above - https://biopython.org/wiki/Alphabet                                                            
                        host_gbk.annotations = {"molecule_type": "DNA"}
                        SeqIO.write(host_gbk, bacteria_filehandle, "genbank")

                        # set the locus tag
                        ppnum += 1
                        pp_gbk = record[phagestart:phageend]
                        pp_gbk.name += f"_PP{ppnum}"
                        # set the Accession
                        pp_gbk.id += f"_PP{ppnum}"
                        # set the definition line
                        pp_gbk.description += f" prophage PP{ppnum} on {contig} from {phagestart} to {phageend}"                                                                      
                        # set the molecule type; required for Biopython v1.78 and above - https://biopython.org/wiki/Alphabet                                                            
                        pp_gbk.annotations = {"molecule_type": "DNA"}
                        SeqIO.write(pp_gbk, phage_filehandle, "genbank")
                        bactstart = phageend + 1

                    host_gbk = record[bactstart:]
                    hostcounter += 1
                    host_gbk.name += f"_region_{hostcounter}"
                    host_gbk.id += f"_region_{hostcounter}"
                    host_gbk.description += f" region {hostcounter} on {contig} from {bactstart} onwards"
                    # set the molecule type; required for Biopython v1.78 and above - https://biopython.org/wiki/Alphabet                                                            
                    host_gbk.annotations = {"molecule_type": "DNA"}
                    SeqIO.write(host_gbk, bacteria_filehandle, "genbank")
                else:
                    SeqIO.write(record, bacteria_filehandle, "genbank")
        
        bacteria_filehandle.close()
        phage_filehandle.close()




