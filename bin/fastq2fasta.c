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
      "fastq2fasta [options] <fastq file> <fasta file>\n"
      "Use  - for either input or output to use STDIN/STDOUT\n"
      "Options:\n"
      "\t-r replace spaces in the fasta header line with '_'\n"
      "\t-n <1/2> add the read number to the header line at the first space\n"
   );
}

int main(int argc, char *argv[]) {

	if ( argc < 2) {
		helpme();
		return 1;
	}
	

	int replace = 0;
	char *n = NULL;

	for (;;) {
		switch(getopt(argc, argv, "rn:")) {
			default:
				helpme();
				return 1;
			case -1:
				break;
			case 'r':
				replace = 1;
				continue;
			case 'n':
				n = optarg;
				continue;
		}
		break;
	}

	if (optind +2 != argc) {
		helpme();
		return 1;
	}

	char* fqf = argv[optind];
	char* faf = argv[optind+1]; 

	FILE *fw;
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


	if (strcmp(faf, "-") != 0)
		fw = fopen(faf, "w");
	else
		fw = stdout;

	seq = kseq_init(fp);
	int l;
	while ((l = kseq_read(seq)) >= 0) {
		if (seq->comment.l == 0) {
			seq->comment.s = malloc(1);
			seq->comment.s = '\0';
		}
		if (replace) {
			for (long unsigned int i = 0; i <= strlen(seq->comment.s); i++) {
				if (seq->comment.s[i] == ' ') {
					seq->comment.s[i] = '_';
				}
			}
			fprintf(fw, ">%s_%s\n%s\n", seq->name.s, seq->comment.s, seq->seq.s);
		}
		else if (n) {
			fprintf(fw, ">%s/%s %s\n%s\n", seq->name.s, n, seq->comment.s, seq->seq.s);
		}
		else {
			fprintf(fw, ">%s %s\n%s\n", seq->name.s, seq->comment.s, seq->seq.s);
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
	fclose(fw);
}


