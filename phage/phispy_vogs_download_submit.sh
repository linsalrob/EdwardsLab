#!/bin/bash
#
# Download the assembly file
#
#      curl -Lo assembly_summary_20220130.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
# Identify which genomes we don't have analysed yet
#
#      mkdir 20220130
#      for BASE in $(grep -v ^# assembly_summary_20220130.txt | cut -f 1); do if [ ! -e phispy/${BASE:0:9}/${BASE:0:13}/VOGS/ ]; then echo $BASE; fi; done > 20220130/needed.all
#
# split those files into separate files with 1000 entries per file
#      mkdir 20220130/needed
#      cd 20220130/needed
#      split -a 4 --numeric-suffixes=1 -l 1000 ../needed.all
#      cd ../../
#
# This makes a few hundred files. In my case 533
# 
#      rm -rf sge_err sge_out; mkdir sge_err sge_out; 
#      qsub -cwd -q smallmem -o sge_out -e sge_err -V -t 1-533:1 ./phispy_vogs_download_submit.sh

NEED=0000$SGE_TASK_ID
NEED=${NEED:(-4)}
snakemake -s phispy_vogs_download.snakefile --config filelist=20220130/needed/x$NEED gbk=20220130/gbk output=20220130/phispy assembly=assembly_summary_20220130.txt  --profile sge
