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
#include <sys/stat.h>
#include <dirent.h>

KSEQ_INIT(gzFile, gzread)

#define table_size 10000

void helpme() {
  fprintf(stderr,
      "count_fastq <fastq file or directory>\n"
      "\tUse  - to read from STDIN\n"
      "\tPass a directory to process all fastq files in that directory\n"
   );
}

// count the sequences in a single file
int count_file(char *fqf) {
	kseq_t *seq;

	// fprintf(stderr, "Reading %s\n", fqf);

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
	unsigned long total = 0;
	unsigned long count = 0;
	while ((l = kseq_read(seq)) >= 0) {
		total += strlen(seq->seq.s);
		count++;
	}
	printf("%s\t%li\t%li\n", fqf, count, total);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

int main(int argc, char *argv[]) {

	if ( argc < 2) {
		helpme();
		return 1;
	}
	

	// test if argv[1] is a directory
	struct stat sb;

	if (stat(argv[1], &sb) == 0 && S_ISDIR(sb.st_mode))
	{
		DIR *d;
		struct dirent *dir;
		d = opendir(argv[1]);
		if (d) {
			int r = 0;
			while ((dir = readdir(d)) != NULL) {
				char filepath[strlen(argv[1]) + strlen(dir->d_name) + 1]; 
				strcpy(filepath, argv[1]);
				strcat(filepath, "/");
				strcat(filepath, dir->d_name);

				struct stat sbf;
				if (stat(filepath, &sbf) == 0 && S_ISREG(sbf.st_mode))
					r += count_file(filepath);
				// else 
				// 	fprintf(stderr, "Skipped %s because it is not a regular file. Its type is %o (last 3 numbers are permissions o, g, w; regular file is 0100000)\n", filepath, sbf.st_mode);

				/*
				 * dir->d_type is often returning 0 (DT_UNKNOWN) for valid files, and I think it is a 
				 * filesystem issue. See https://stackoverflow.com/a/47079705
				 if (dir->d_type == DT_REG) {
					char filepath[strlen(argv[1]) + strlen(dir->d_name) + 1]; 
					strcpy(filepath, argv[1]);
					strcat(filepath, "/");
					strcat(filepath, dir->d_name);
					r += count_file(filepath);
				} else 
					fprintf(stderr, "Skipped %s because it is not a regular file. Its type is %d\n", dir->d_name, dir->d_type);
				*/
			}
			closedir(d);
			return r;
		}
	}

	return count_file(argv[1]);

}


