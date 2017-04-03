#!/bin/bash

#############################################
#
#  Compare all the reads in the metagneome
#
##############################################

IFS='
'

################# change these file names 
#
# this is the output location
ODIR=/vol1/bam.crassphage
# this is the name of the bowtie index, without the suffixes
BOWTIEIDX=crassphage
#
################ you don't need to change anything below this



if [ ! -n "$1" ]
then
        echo "Usage: `basename $0` <file with list of SRA IDs>";
        exit $E_BADARGS
fi


touch START$$
DIR=fastq$$

if [ ! -e $ODIR ]; then mkdir $ODIR; fi

for i in $(cat $1); do 
	if [ ! -e $ODIR/$i.bam ]; then
		echo $i;
		fastq-dump --outdir $DIR --skip-technical  --readids --read-filter pass --dumpbase --split-files --clip $i
		READS=$(ls $DIR/* | tr \\n \, | sed -e 's/,$//')
		bowtie2 -p 6 -q --no-unal -x $BOWTIEIDX -U $READS | samtools view -bS - | samtools sort - $ODIR/$i
		samtools index $ODIR/$i.bam
		rm -rf $DIR
		rm -f $HOME/ncbi/public/sra/$i.sra*
	fi
done
