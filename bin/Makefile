# A simple Makefile to compile the C code here
#

FLAGS := $(FLAGS) -Wall
override CFLAGS += $(shell pkg-config --cflags --libs python3-embed)

all:
	gcc $(FLAGS) -o fastq_avqual fastq_avqual.c -lz
	gcc $(FLAGS) -o count_fasta countfasta.c -lz
	c++ $(FLAGS) -o fastq2fasta fastq2fasta.cpp 


clean:
	rm -f fastq_avqual count_fasta fastq2fasta

