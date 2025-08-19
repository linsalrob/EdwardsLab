/* fastq_minmax.c
 * Find the shortest and longest sequences in a FASTQ using kseq.h
 * Compile: cc -O2 -o fastq_longest_shortest fastq_longest_shortest.c -lz
 * Written by ChatGPT 5
 * 19/8/2025
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <zlib.h>
#include "kseq.h"           
KSEQ_INIT(gzFile, gzread)

/* ANSI C doesn't guarantee strdup; provide a tiny replacement. */
static char *xstrdup(const char *s) {
    size_t n;
    char *p;
    if (s == NULL) return NULL;
    n = strlen(s) + 1;
    p = (char *)malloc(n);
    if (p == NULL) return NULL;
    memcpy(p, s, n);
    return p;
}

static void die(const char *msg) {
    fprintf(stderr, "Error: %s\n", msg);
    exit(EXIT_FAILURE);
}

int main(int argc, char *argv[]) {
    gzFile fp = NULL;
    kseq_t *seq = NULL;
    long nreads = 0;

    /* Trackers for min and max */
    unsigned long min_len = 0, max_len = 0;
    char *min_name = NULL, *max_name = NULL;
    char *min_seq = NULL,  *max_seq = NULL;

    if (argc != 2) {
        fprintf(stderr, "Usage: %s <fastq.gz|fastq|- for stdin>\n", argv[0]);
        return EXIT_FAILURE;
    }

    if (strcmp(argv[1], "-") == 0) {
        /* Read from stdin */
        fp = gzdopen(fileno(stdin), "r");
        if (fp == NULL) die("failed to open stdin with gzdopen()");
    } else {
        fp = gzopen(argv[1], "r");
        if (fp == NULL) {
            fprintf(stderr, "Error: cannot open '%s': %s\n", argv[1], strerror(errno));
            return EXIT_FAILURE;
        }
    }

    seq = kseq_init(fp);
    if (seq == NULL) die("kseq_init failed");

    while (kseq_read(seq) >= 0) {
        unsigned long L = (unsigned long)seq->seq.l; /* cast for ANSI printf */
        ++nreads;

        if (nreads == 1) {
            /* First record initializes both min and max */
            min_len = max_len = L;
            min_name = xstrdup(seq->name.s);
            max_name = xstrdup(seq->name.s);
            min_seq  = xstrdup(seq->seq.s);
            max_seq  = xstrdup(seq->seq.s);
            if (!min_name || !max_name || !min_seq || !max_seq) die("out of memory");
        } else {
            if (L < min_len) {
                char *new_name = xstrdup(seq->name.s);
                char *new_seq  = xstrdup(seq->seq.s);
                if (!new_name || !new_seq) die("out of memory");
                free(min_name); free(min_seq);
                min_name = new_name; min_seq = new_seq; min_len = L;
            }
            if (L > max_len) {
                char *new_name = xstrdup(seq->name.s);
                char *new_seq  = xstrdup(seq->seq.s);
                if (!new_name || !new_seq) die("out of memory");
                free(max_name); free(max_seq);
                max_name = new_name; max_seq = new_seq; max_len = L;
            }
        }
    }

    if (nreads == 0) {
        fprintf(stderr, "No reads found.\n");
        kseq_destroy(seq);
        gzclose(fp);
        return EXIT_FAILURE;
    }

    /* Report */
    printf("File\tReads Processed\tShortest Read\tLongest Read\n");
    printf("%s\t%ld\t%lu (%s)\t%lu (%s)\n", argv[1], nreads, min_len, min_name ? min_name : "(unknown)", max_len, max_name ? max_name : "(unknown)");

    /*
    printf("Total reads processed: %ld\n\n", nreads);

    printf("Shortest read:\n");
    printf("  ID: %s\n", min_name ? min_name : "(unknown)");
    printf("  Length: %lu\n", min_len);
    printf("  Sequence:\n%s\n\n", min_seq ? min_seq : "");

    printf("Longest read:\n");
    printf("  ID: %s\n", max_name ? max_name : "(unknown)");
    printf("  Length: %lu\n", max_len);
    printf("  Sequence:\n%s\n", max_seq ? max_seq : "");
    */

    /* Cleanup */
    free(min_name); free(min_seq);
    free(max_name); free(max_seq);
    kseq_destroy(seq);
    gzclose(fp);

    return EXIT_SUCCESS;
}

