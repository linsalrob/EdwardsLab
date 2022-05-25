#!/bin/bash

DATE=20220525
ASS=$DATE/assembly_summary_$DATE.txt.gz
VOGS=/home3/redwards/VOGs/VOGs.hmm 

NEED=0000$SGE_TASK_ID
NEED=${NEED:(-4)}

snakemake -s ~/GitHubs/EdwardsLab/phage/phispy_vogs_download.snakefile --config filelist=$DATE/needed/x$NEED gbk=$DATE/gbk output=$DATE/phispy assembly=$ASS vogs=$VOGS --profile sge 
