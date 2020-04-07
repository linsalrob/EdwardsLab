/*
 * Print some sequences out of a fastq file
 *
 * This is written for search SRA so I can easily limit a fastq file
 * to 10,000 reads. I could probably do that with samtools, but this
 * is lightweight, and I needed practice!
 *
 * Created by Rob on 4/6/20.
 * Leans heavily on the kseq library from Heng Li: https://github.com/lh3/readfq
 *
 * If you have the kseq.h library installed, this should compile with 
 * no warnings (see below) with:
 *
 * gcc -o print_fastq print_fastq.c -Wall -lz
 *
 * Usage:
 *     print_fastq -n <number of sequences> <fastq file>
 *
 * With the -v option you get the number of reads written on stderr. You can redirect
 * the stdout to gzip and stderr to a file like so:
 *     print_fastq -n 5 -v test.fastq.gz 2>> num_reads.tsv | gzip - > test.out.gz
 *
 * If you do not have the kseq.h library you can download it from
 * https://github.com/lh3/readfq and then compile (eg. if it is in ../lib)
 * gcc -o print_fastq print_fastq.c -Wall -I../lib -lz
 *
 * For some versions of gcc you may need to add
 * -std=c99  if you get error: ‘for’ loop initial declarations are only allowed in C99 mode
 * -D_XOPEN_SOURCE=700  if you get warning: implicit declaration of function ‘fileno’
 *
 *  to give:
 *  gcc -o print_fastq print_fastq.c -Wall -I../lib -lz -std=c99 -D_XOPEN_SOURCE=700
 *
 */

#include <zlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("%s [-n number of reads] [-v] [fastq file]\n", argv[0]);
        return 1;
    }

    char *filename;
    int verbose = 0;
    int nseqs = 0;
    struct stat buffer;
    for (int i=0; i<argc; i++) {
        if (strcmp(argv[i], "-n") == 0) {
            if (sscanf(argv[i + 1], "%i", &nseqs) != 1) {
                fprintf(stderr, "Error: -n must be followed by a number of sequences to print");
                return 1;
            }
        }
        else if (strcmp(argv[i], "-v") == 0)
            verbose = 1;
        else if (stat(argv[i],&buffer) == 0)
            filename = argv[i];
    }

    gzFile fp;
    kseq_t *seq;
    FILE * file = fopen(filename, "rb");
    fp = gzdopen(fileno(file), "r");
    seq = kseq_init(fp);
    int counter = 0;
    while (kseq_read(seq) >= 0) {
        printf("@%s %s\n", seq->name.s, seq->comment.s);
        printf("%s\n+\n%s\n", seq->seq.s, seq->qual.s);
        counter++;
        if (counter == nseqs) {
            if (verbose) { fprintf(stderr, "%s\t%i\n", filename, counter); }
            return 0;
        }
    }
    kseq_destroy(seq);
    gzclose(fp);
    if (verbose) {fprintf(stderr, "%s\t%i\n", filename, counter);}
    return 0;
}
