/*
 *
 * Change the character that separates the pair read number at the end of the 
 * read name.
 *
 * Some fastq files from SRA use .1 & .2 at the end of the read name (before the first
 * space) to indicate read pairs instead of /1 & /2. This changes that
 *
 * Created by Rob on 4/6/20.
 * Leans heavily on the kseq library from Heng Li: https://github.com/lh3/readfq
 *
 * If you have the kseq.h library installed, this should compile with 
 * no warnings (see below) with:
 *
 * gcc -o change_fastq_pair_symbol change_fastq_pair_symbol.c -Wall -lz
 *
 * Usage:
 *     change_fastq_pair_symbol  <fastq file>
 *
 * If you do not have the kseq.h library you can download it from
 * https://github.com/lh3/readfq and then compile (eg. if it is in ../lib)
 *
 * For some versions of gcc you may need to add
 * -std=c99  if you get error: ‘for’ loop initial declarations are only allowed in C99 mode
 * -D_XOPEN_SOURCE=700  if you get warning: implicit declaration of function ‘fileno’
 *
 *  to give:
 *  gcc -o change_fastq_pair_symbol change_fastq_pair_symbol.c -Wall -I../lib -lz -std=c99 -D_XOPEN_SOURCE=700
 *
 */

#include <zlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[]) {
    if (argc < 1) {
        printf("%s [-v] [fastq file]\n", argv[0]);
        return 1;
    }

    char *filename;
    int verbose = 0;
    int nseqs = 0;
    struct stat buffer;
    for (int i=0; i<argc; i++) {
        if (strcmp(argv[i], "-v") == 0)
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
	// get the last two characters of the name
	int l = strlen(seq->name.s);
	
	seq->name.s[l-2] = '/';

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
