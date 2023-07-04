#!/bin/bash

###############################################################
#                                                             #
# submit a directory of R1 and R2 fasta files for mmseqs easy taxonomy
#                                                             #
#                                                             #
###############################################################



if [[ $# < 3 ]]; then echo `basename $0` "<database UniRef50|GTDB> <output drectory for all files> <directory  of fasta files>" >&2;
	exit 2;
fi

DB=$1

if [[ $DB == "NR" ]]; then
	echo "You need at least 1.5TB ram for the NR comparison. Please used the dedicated slurm script" >&2;
	exit 3;
fi

if [[ $DB != "UniRef50" && $DB != "GTDB" ]]; then
	echo "Please use one of UniRef50, NR, or GTDB for the current database" >&2; 
	exit 2
fi

OUTDIR=$2
mkdir --parents $OUTDIR

FASTA=$3

for R1 in $(find $FASTA -name \*R1\* -printf "%f\n"); do
	R2=${R1/R1/R2}
	OUT=$(echo $R1 | sed -e 's/.R1.*$//')
	echo "submitting: sbatch /home/edwa0468/GitHubs/EdwardsLab/process_JCJ_metagenomes/mmseqs_easy_taxonomy.slurm $DB $OUTDIR/$OUT $FASTA/$R1 $FASTA/$R2";
	sbatch /home/edwa0468/GitHubs/EdwardsLab/process_JCJ_metagenomes/mmseqs_easy_taxonomy.slurm $DB $OUTDIR/$OUT $FASTA/$R1 $FASTA/$R2;
done
