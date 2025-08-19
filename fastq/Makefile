# A simple Makefile to compile the C code here
#

FLAGS := $(FLAGS) -Wall -O3
override CFLAGS += $(shell pkg-config --cflags --libs python3-embed)

all:
	gcc $(FLAGS) -o fastq_longest_shortest fastq_longest_shortest.c -lz

clean:
	rm -f fastq_longest_shortest

