#!/bin/bash
#SBATCH --job-name=mmseqs_submit
#SBATCH --time=0-12
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -o mmseqs_submit-%j.out
#SBATCH -e mmseqs_submit-%j.err
#SBATCH --partition=short

###############################################################
#                                                             #
# submit a directory of R1 and R2 fasta files for mmseqs easy taxonomy
# this is so we can submit it as a slurm job and use dependency
# sbatch --dependency=afterok:$JOB job.slurm
#                                                             #
#                                                             #
###############################################################

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FAFILEEND, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh

FASTA=fasta

if [[ ! -e $FASTA ]]; then
	echo "ERROR: Directory $FASTA not found" >&2;
	exit 2;
fi

DB=UniRef50

OUTDIR=mmseqs
mkdir --parents $OUTDIR


for R1 in $(find $FASTA -name \*R1\* -printf "%f\n"); do
	R2=${R1/R1/R2}
	OUT=${R1/$FAFILEEND/}
	if [[ ! -e $OUTDIR/$OUT ]]; then
		echo "submitting: sbatch /home/edwa0468/GitHubs/EdwardsLab/process_EK_metagenomes/mmseqs_easy_taxonomy.slurm $DB $OUTDIR/$OUT $FASTA/$R1 $FASTA/$R2";
		sbatch /home/edwa0468/GitHubs/EdwardsLab/process_EK_metagenomes/mmseqs_easy_taxonomy.slurm $DB $OUTDIR/$OUT $FASTA/$R1 $FASTA/$R2;
	else 
		echo "Output directory $OUTDIR/$OUT already exists, so did not submit" >&2;
	fi
done
