#
#  fastq2crassphage_anthill.sh
#
# Compare all the fastq reads to crAssphage using bowtie and index files with significant read alignments.
#
# To submit this to the cluster, you need fastq_paired_files.txt that is created by fastq_pairs.py
# then use wc to count the number of lines in that file, and qsub to submit it to the cluster:
#
#            wc -l fastq_paired_files.txt
# 182787
#            qsub -cwd -o cpsge/ -e cpsge/ -t 1-75000:1 ./fastq2crassphage_anthill.sh
#            qsub -cwd -o cpsge/ -e cpsge/ -t 75001-150000:1 ./fastq2crassphage_anthill.sh 
#            qsub -cwd -o cpsge/ -e cpsge/ -t 150001-182787:1 ./fastq2crassphage_anthill.sh 
#
#
# Rob Edwards, Dec 2015
#
# updated to (1) handles gzipped files and (2) handles paired end reads


# This is the location of our indexed crAssphage sequence.
export BOWTIE2_INDEXES=/home3/redwards/phage/crAssphage/

# read a line 
FQ=$(head -n $SGE_TASK_ID fastq_paired_files.txt | tail -n 1)

# echo  the ID to stderr, so that we know what the fastq file was if there is any bowtie output to stderr
>&2 echo $FQ


# the revised fastq_pairs.py prints out the complete line we need so we don't need to split on spaces any more
# this allows it to correct for uneven lengthed fastq files

# split this on spaces. The second option will be a filename
IFS=' ' read -r -a ARRAY <<< $FQ

OUT=$(echo ${ARRAY[1]} | sed -e 's/.clean//; s/.unique//; s/.fastq.gz//; s/.fastq$//; s/fastq/crassphage/;  s/_[0-9]\+$//')
>&2 echo "OUTPUT TO $OUT.sam"

if [ ! -e $OUT.sam ]; then
	$HOME/bin/bowtie/current/bowtie2 -q -x JQ995537 --no-unal $FQ -S $OUT.sam
else
	echo "Not overwriting $OUT.sam";
fi;

# If the file is over 1000 bytes it has something interesting. Convert it to a bam file and index it
FILESIZE=$(wc -c <"$OUT.sam")
if [ "$FILESIZE" -gt "1000" ]; then 
 	$HOME/bin/samtools/samtools-1.2/samtools view -bS $OUT.sam | $HOME/bin/samtools/samtools-1.2/samtools sort - $OUT;
 	$HOME/bin/samtools/samtools-1.2/samtools index $OUT.bam;
fi

