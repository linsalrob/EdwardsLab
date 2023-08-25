#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

void helpme() {
  fprintf(stderr,
      "fasta_split [options] <interleaved fasta file> <fasta R1 file> <fasta R2 file>\n"
   );
}

int main(int argc, char *argv[])
{
	if ( argc < 3) {
		helpme();
		return 1;
	}
	
	int c;
	gzFile fp;
	kseq_t *ks;

	while ((c = getopt(argc, argv, "")) >= 0);

	if (optind == argc) {
		helpme();
		return 1;
	}

	fp = gzopen(argv[optind], "r");
	if (fp == NULL) {
		fprintf(stderr, "ERROR: There was a problem opening %s for reading\n", argv[optind]);
		return 1;
	}

	ks = kseq_init(fp);

	FILE *fr1 = fopen(argv[optind+1], "w");
	FILE *fr2 = fopen(argv[optind+2], "w");

	if (fr1 == NULL || fr2 == NULL) {
		fprintf(stderr, "ERROR: There was a problem opening %s or %s for writing\n", argv[optind+1], argv[optind+2]);
		return 1;
	}

	while (kseq_read(ks) >= 0) {
		if  (ks->name.s[ks->name.l - 1] == '1')
			fprintf(fr1, ">%s %s\n%s\n",  ks->name.s, ks->comment.s, ks->seq.s);
		else
			fprintf(fr2, ">%s %s\n%s\n",  ks->name.s, ks->comment.s, ks->seq.s);
	}

	kseq_destroy(ks);
	gzclose(fp);
	
	fclose(fr1);
	fclose(fr2);

	return 0;
}
