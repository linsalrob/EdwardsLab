/* Average the quality scores in a fastq file 
 *
 * Rob Edwards, 10/12/21
 *
 */


#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

#define table_size 10000

void helpme() {
  fprintf(stderr,
      "fastq2fasta [options] <fastq file>\n"
      "Use  - to read from STDIN\n"
   );
}

int main(int argc, char *argv[]) {

	if ( argc < 2) {
		helpme();
		return 1;
	}
	

	char* fqf = argv[1];

	kseq_t *seq;


	FILE *instream = NULL;
 
	if (strcmp(fqf, "-") !=0) {
		instream = fopen(fqf, "r");
		if (instream == NULL) {
			fprintf(stderr, "ERROR: Can't open %s\n", fqf);
			return 1;
		}
	}
	else
		instream = stdin;

	gzFile fp = gzdopen(fileno(instream), "r");
	if (fp == NULL) {
		fprintf(stderr, "ERROR: There was a problem opening the stream to %s\n", fqf);
		return 1;
	}



	seq = kseq_init(fp);
	int l;
	int total = 0;
	int count = 0;
	while ((l = kseq_read(seq)) >= 0) {
		total += strlen(seq->seq.s);
		count++;
	}
	printf("%s\t%i\t%i\n", fqf, count, total);
	kseq_destroy(seq);
	gzclose(fp);
}


