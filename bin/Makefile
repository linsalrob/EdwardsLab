# A simple Makefile to compile the C code here
#

FLAGS := $(FLAGS) -Wall
override CFLAGS += $(shell pkg-config --cflags --libs python3-embed)

all:
	gcc $(FLAGS) -o fastq_avqual fastq_avqual.c -lz
	gcc $(FLAGS) -o count_fasta count_fasta.c -lz
	# c++ $(FLAGS) -o fq2fa fastq2fasta.cpp 
	gcc $(FLAGS) -o fastq2fasta fastq2fasta.c -lz
	gcc $(FLAGS) -o count_fastq count_fastq.c -lz

clean:
	rm -f fastq_avqual count_fasta fastq2fasta count_fastq

