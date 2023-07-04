#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

/*
 * fastg2gfa downloaded from GitHub 2/7/2023
 * https://raw.githubusercontent.com/lh3/gfa1/master/misc/fastg2gfa.c
 *
 * Please use the original and don't use this version.
 * Please cite https://github.com/lh3/gfa1
 */

int main(int argc, char *argv[])
{
	int c;
	gzFile fp;
	kseq_t *ks;

	while ((c = getopt(argc, argv, "")) >= 0);

	if (optind == argc) {
		fprintf(stderr, "Usage: fastg2gfa <in.fastg>\n");
		return 1;
	}

	fp = gzopen(argv[optind], "r");
	assert(fp);
	ks = kseq_init(fp);

	while (kseq_read(ks) >= 0) {
		char *p, *s = ks->name.s;
		int is_comp;
		for (p = s; *p && *p != ':' && *p != ';'; ++p);
		c = *p, *p = 0;
		is_comp = (p > s && *(p-1) == '\''); // if we are looking at a complement segment
		if (is_comp) *(p-1) = 0;
		if (!is_comp)
			printf("S\t%s\t%s\tLN:i:%ld\n", s, ks->seq.s, (long)strlen(ks->seq.s));
		if (c == ':') { // have neighbors
			char *q = p + 1;
			do {
				int is_comp2 = 0;
				for (p = q; *p && *p != ',' && *p != ';'; ++p);
				c = *p, *p = 0;
				is_comp2 = (p > q && *(p-1) == '\'');
				if (is_comp2) *(p-1) = 0;
				printf("L\t%s\t%c\t%s\t%c\t0M\n", s, "+-"[!!is_comp], q, "+-"[!!is_comp2]);
				q = p + 1;
			} while (c != 0 && c != ';');
		}
	}

	kseq_destroy(ks);
	gzclose(fp);
	return 0;
}
