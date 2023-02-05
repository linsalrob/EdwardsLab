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
		printf("Usage: %s <fastq file (use - to read from STDIN)>\n", argv[0]);
		return 1;
	}
	gzFile fp;
	kseq_t *seq;
	int c=0;
	long total=0;
	long n=0;

	fp = gzopen(argv[1], "r");
	seq = kseq_init(fp);
	int l;
	while ((l = kseq_read(seq)) >= 0) {
		c++;
		for (int i = 0; i < strlen(seq->qual.s); i++) {
			total+= (int) seq->qual.s[i];
			n++;
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
	printf("File\tNumber of sequences\tTotal bp\tTotal quality\tAverage quality\n");
	printf("%s\t%d\t%ld\t%ld\t%ld\n", argv[1], c, n, total, total/n);

}

