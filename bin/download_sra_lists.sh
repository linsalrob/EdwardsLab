#!/bin/bash

USAGE="
$0 [previous downloaded files (optional)] [file with list of SRA IDs] [destination directory]

Download a set of SRA IDs into a destination directory

You may also provide another directory that has the files you have already downloaded
and we will skip any duplicates. Duplicates will be detected by finding the first part
of the name up until a '_' or '.' so SRR061207_2.fastq.gz or SRR061207.fastq will both
be considered as SRR061207.
"

E_BADARGS=65
E_NOFILE=66

DEBUG=0

if [ "$#" -eq 3 ]; then
	PREVDIR=$1
	INPUT=$2
	OUTDIR=$3
elif [ "$#" -eq 2 ]; then
	PREVDIR=""
	INPUT=$1
	OUTDIR=$2
else
	echo "$USAGE"
	exit $E_BADARGS
fi

if [ $DEBUG -gt 0 ]; then
	echo "Previous files are in $PREVDIR"
	echo "Input text file with list of SRAs is in $INPUT"
	echo "Output will be to $OUTDIR"
fi

if [ ! -e $INPUT ]; then
	echo "$INPUT does not exist"
	exit $E_NOFILE
fi

if [ ! -e $OUTDIR ]; then
	mkdir -p $OUTDIR
fi

declare -A ALREADYHAVE
if [ -d $PREVDIR ]; then
	echo "Reading files from $PREVDIR"
	
	# this is done using an associative array, of course
	for SRR in $(find $PREVDIR -type f -printf "%f\n"  | sed -E 's/[_\.].*$//'); do ALREADYHAVE[$SRR]=1; done
	
	if [ $DEBUG -gt 0 ]; then
		# print all previously found entries
		# for k in "${!ALREADYHAVE[@]}"; do echo "PREVIOUS METAGENOME: $k --> ${ALREADYHAVE[$k]}"; done
		for k in "${!ALREADYHAVE[@]}"; do echo "PREVIOUS METAGENOME: $k"; done
	else
		echo "There are ${#ALREADYHAVE[@]} previous metagenomes that we have remembered to ignore!"
	fi
fi

for SRR in $(cat $INPUT); do
	if [[ ${ALREADYHAVE[$SRR]} == 1 ]]; then
		if [ $DEBUG -gt 0 ]; then
			echo "Skipping $SRR";
		fi
	else
		fastq-dump --outdir $OUTDIR --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip -N 1000 -X 101000 $SRR
	fi
done

