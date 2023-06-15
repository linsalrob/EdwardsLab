#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "kseq.h"

KSEQ_INIT(gzFile, gzread);

// compile with:
// gcc -I../include -o filter_fastq_length ./filter_fastq_length.c -lz


int main(int argc, char* argv[]) {
	if (argc < 2) {
		fprintf(stderr,  "%s <fastq file> <minimum length>\n", argv[0]);
		exit(1);
	}

	gzFile fp;
	kseq_t *seq;
	fp = gzopen(argv[1], "r");
	seq = kseq_init(fp);
	int l;
	int minlen = atoi(argv[2]);
	fprintf(stderr, "Filtering %s. Reads shorter than %d will be ignored\n", argv[1], minlen);
	while ((l = kseq_read(seq)) >= 0) {
		if (seq->seq.l > minlen)
			printf("@%s %s\n%s\n+\n%s\n", seq->name.s, seq->comment.s, seq->seq.s, seq->qual.s);
	}
	return 0;
}

