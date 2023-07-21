#!/bin/bash

###############################################################
#                                                             #
# mmseqs easy taxonomy against UniRef50                       #
#                                                             #
# Note that you need a fasta file                             #
# Also make sure that the reads are labeled /1 /2             #
# (without a space)                                           #
#                                                             #
# sbatch mmseqs_easy_taxonomy.slurm                           #
#                                                             #
#                                                             #
###############################################################


#SBATCH --job-name=mmseqsLT
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH -o mmseqsLTE-%j.out
#SBATCH -e mmseqsLTE-%j.err

eval "$(conda shell.bash hook)"
conda activate bioinformatics
set -euo pipefail



# Output file is the base name. There will be several files with that name, so we write to $BGFS and then copy everything to the directory
# with that name once done

echo "Start " `date` >&2;

if [[ $# < 3 ]]; then echo `basename $0` "<database Uniref50|GTDB> <output file> <fasta files>" >&2;
	exit 2;
fi

DB=$1
shift

if [[ $DB == "NR" ]]; then
	echo "You need at least 1.5TB ram for the NR comparison. Please used the dedicated slurm script" >&2;
	exit 3;
fi

if [[ $DB != "UniRef50" && $DB != "GTDB" ]]; then
	echo "Please use one of UniRef50, NR, or GTDB for the current database" >&2; 
	exit 2
fi

OUTPUT=$1
shift

TMPOUTPUT=`basename $OUTPUT`
mkdir $BGFS/output
echo "Temporary output is in $BGFS/output/$TMPOUTPUT*"

FASTA="$@"


SOURCE=/scratch/user/edwa0468/Databases/mmseqs/$DB

if [[ ! -e $SOURCE ]]; then 
	echo "FATAL: $SOURCE does not exist. Please check and try again";
	exit 2;
fi



mkdir -p $BGFS/$DB
rsync -a $SOURCE/ $BGFS/$DB

echo "Copied " `date` >&2;

echo "Running mmseqs easy-taxonomy $FASTA $BGFS/$DB/$DB $BGFS/output/$TMPOUTPUT tmp --start-sens 1 --sens-steps 3 -s 7 --threads 32" 2>&1;

mmseqs easy-taxonomy $FASTA $BGFS/$DB/$DB $BGFS/output/$TMPOUTPUT $(mktemp -d -p $BGFS) --start-sens 1 --sens-steps 3 -s 7 --threads 32

find $BGFS/output -type f  | parallel -j 32 gzip

echo "Copying results from $BGFS/output* to $OUTPUT";

mkdir -p $OUTPUT
rsync -av $BGFS/output/ $OUTPUT/

echo "Fin " `date` >&2;

