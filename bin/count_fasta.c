/* 
 * count the sequences in a fasta file
 *
 * just prints the number of sequences and number of bases at the moment
 */

#include <stdlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "zlib.h"
#include <dirent.h> 
#include <getopt.h>

#if defined(WIN32)
#  define DIR_SEPARATOR '\\'
#else
#  define DIR_SEPARATOR '/'
#endif


void helpme() {
	fprintf(stderr,
			"countfasta [-p print all read lengths] [-f fasta file] [-d directory]\n"
	       );
}


void join_path(char *destination, const char *path1, const char *path2) {
 // this comes from https://stackoverflow.com/questions/3142365/combine-directory-and-file-path-c
  if (path1 && *path1) {
    int len = strlen(path1);
    strcpy(destination, path1);

    if (destination[len - 1] == DIR_SEPARATOR) {
      if (path2 && *path2) {
        strcpy(destination + len, (*path2 == DIR_SEPARATOR) ? (path2 + 1) : path2);
      }
    }
    else {
      if (path2 && *path2) {
        if (*path2 == DIR_SEPARATOR)
          strcpy(destination + len, path2);
        else {
          destination[len] = DIR_SEPARATOR;
          strcpy(destination + len + 1, path2);
        }
      }
    }
  }
  else if (path2 && *path2)
    strcpy(destination, path2);
  else
    destination[0] = '\0';
}

void count_file(char *filename, int print_all, int verbose) {
	gzFile fp;

	int MAXLINELEN = 100000;
	
	char line[MAXLINELEN];

	if ((fp = gzopen(filename, "r")) == NULL) {
		if (verbose)
			fprintf(stderr, "Can't open file %s\n", filename);
		exit(1);
	}

	size_t len = 0;
	int n = 0;
	char *lastid = NULL;
	int this_seq_len = 0;
	while ((gzgets(fp, line, MAXLINELEN)) != NULL) {
		if ((int) line[0] == 62) { // not sure why I'm using an ascii comparison, but I'm thinking ascii at the moment
			n++;
			if (lastid && print_all) {
				lastid[strcspn(lastid, "\r\n")] = 0;
				printf("%s\t%i\n", lastid, this_seq_len);
			}
			lastid = line;
			this_seq_len = 0;
		}
		else {
			len += strlen(line);
			this_seq_len += strlen(line);
		}
	}
	if (print_all)
		printf("%s\t%i\n", lastid, this_seq_len);
	printf("%s\t%i\t%zu\n", filename, n, len);
}

void read_directory(char *dirname, int print_all, int verbose) {
	DIR *d;
	struct dirent *dir;
	d = opendir(dirname);
	if (d) {
		while ((dir = readdir(d)) != NULL) {
			if (strcmp(dir->d_name, ".") == 0 || strcmp(dir->d_name, "..") == 0)
				continue;
			char full_path[strlen(dirname) + strlen(dir->d_name) + 2];
			join_path(full_path, dirname, dir->d_name);
			count_file(full_path, print_all, verbose);
		}
		closedir(d);
	}
}


int main(int argc, char *argv[]) {
	int verbose = 0;
	int done = 0;
	int print_all = 0;
	for (;;) {
		switch(getopt(argc, argv, "pd:f:v:")) {
			default:
				helpme();
				return 1;
			case -1:
				break;
			case 'v':
				verbose = 1;
			case 'd':
				read_directory(optarg, print_all, verbose);
				done = 1;
				continue;
			case 'p':
				print_all = 1;
				continue;
			case 'f':
				count_file(optarg, print_all, verbose);
				done = 1;
				continue;
		}
		break;
	}
	if (done == 0) {
		helpme();
		return 1;
	}
}


