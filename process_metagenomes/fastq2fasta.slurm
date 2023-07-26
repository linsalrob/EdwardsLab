#!/bin/bash
#SBATCH --job-name=fastq2fasta
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -o faing-%j.out
#SBATCH -e faing-%j.err

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate bioinformatics


if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

# Note: according to Scott, writing to BGFS and then copying to /scratch is going to be a lot faster
# than just writing to /scratch.

# this requires fastq2fasta.c be compiled which has built in reading gzip support

extract_R1 () {
	FQ=$1;
	FB=`basename $FQ .gz`;
	FA=${FB/fastq/fasta}
	fastq2fasta -n 1 $FQ - | gzip -c > $BGFS/$FA.gz;
	cp -t fasta $BGFS/$FA.gz
}

extract_R2 () {
	FQ=$1;
	FB=`basename $FQ .gz`;
	FA=${FB/fastq/fasta}
	fastq2fasta -n 2 $FQ - | gzip -c > $BGFS/$FA.gz;
	cp -t fasta $BGFS/$FA.gz
}

export -f extract_R1
export -f extract_R2
mkdir -p fasta

find -L $HOSTREMOVED -type f -name \*R1\*fastq.gz | parallel -j 16 extract_R1
find -L $HOSTREMOVED -type f -name \*R2\*fastq.gz | parallel -j 16 extract_R2

