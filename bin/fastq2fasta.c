/* Average the quality scores in a fastq file 
 *
 * Rob Edwards, 10/12/21
 *
 */


#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

#define table_size 10000

int main(int argc, char *argv[]) {

	if ( argc < 2) {
		printf("Usage: %s <fastq file (use - to read from STDIN)> <fasta file (use - to write to STDOUT)>\n", argv[0]);
		return 1;
	}
	
	FILE *fw;
	kseq_t *seq;


	FILE *instream = NULL;
 
	if (strcmp(argv[1], "-") !=0) {
		instream = fopen(argv[1], "r");
		if (instream == NULL) {
			fprintf(stderr, "ERROR: Can't open %s\n", argv[1]);
			return 1;
		}
	}
	else
		instream = stdin;
	gzFile fp = gzdopen(fileno(instream), "r");
	if (fp == NULL) {
		fprintf(stderr, "ERROR: There was a problem opening the stream to %s\n", argv[1]);
		return 1;
	}


	if (strcmp(argv[2], "-") != 0)
		fw = fopen(argv[2],"w");
	else
		fw = stdout;

	seq = kseq_init(fp);
	int l;
	while ((l = kseq_read(seq)) >= 0) {
		fprintf(fw, ">%s\n%s\n", seq->name.s, seq->seq.s);
	}
	kseq_destroy(seq);
	gzclose(fp);
	fclose(fw);
}


