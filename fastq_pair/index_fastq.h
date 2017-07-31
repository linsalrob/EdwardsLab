//
// Created by redwards on 7/29/17.
//

#ifndef CEEQLIB_INDEX_FASTQ_H
#define CEEQLIB_INDEX_FASTQ_H

#include <stdbool.h>

struct idloc {
    long int pos;
    char *id;
    bool printed;
};

struct filehash {
    struct filehash *next;
    struct idloc *id;
};

// how long should our lines be. This is a 64k buffer
#define MAXLINELEN 65536
// what is our hash table size. This should probably approach the size
// of our fastq file
#define HASHSIZE 2000000


/*
 * calculate the hash for a fastq sequence
 *
 * This is a simple hash but widely used!
 *
 * we use an unsigned here so that the answer is > 0
 */

unsigned hash (char *s);

int index_file(char *f, char *g);

#endif //CEEQLIB_INDEX_FASTQ_H
