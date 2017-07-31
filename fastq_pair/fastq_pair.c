/*
 * Pair a fastq file using a quick index
 */

#include "fastq_pair.h"
#include "robstr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// our hash for the fastq ids
struct filehash *ids[HASHSIZE] = {NULL};

int DEBUG = 1;

int index_file(char *left_fn, char *right_fn) {
    FILE *lfp;

    char *line = malloc(sizeof(char) * MAXLINELEN + 1);

    if ((lfp = fopen(left_fn, "r")) == NULL) {
        fprintf(stderr, "Can't open file %s\n", left_fn);
        exit(1);
    }

    char *aline; /* this variable is not used, it suppresses a compiler warning */

    long int nextposition = 0;

    /*
     * Read the first file and make an index of that file.
     */
    while ((aline = fgets(line, MAXLINELEN, lfp)) != NULL) {
        struct idloc *newid;
        newid = (struct idloc *) malloc(sizeof(*newid));

        if (newid == NULL) {
            fprintf(stderr, "Can't allocate memory for new ID pointer\n");
            return 0;
        }

        line[strcspn(line, " ")] = '\0';


        /*
         * Figure out what the match mechanism is. We have three examples so
         *     i.   using /1 and /2
         *     ii.  using /f and /r
         *     iii. just having the whole name
         *
         * If there is a /1 or /2 in the file name, we set that part to null so the string is only up
         * to before the / and use that to store the location.
         */

        char lastchar = line[strlen(line)-1];
        if ('1' == lastchar || '2' == lastchar || 'f' == lastchar ||  'r' == lastchar)
            line[strlen(line)-1] = '\0';


        newid->id = dupstr(line);
        newid->pos = nextposition;
        newid->printed = false;

        struct filehash *ptr;
        ptr = (struct filehash *) malloc(sizeof(*ptr));
        if (ptr == NULL) {
            fprintf(stderr, "Can't allocate memory for pointer\n");
            return -1;
        }

        ptr->id = newid;

        unsigned hashval = hash(newid->id);
        ptr->next = ids[hashval];
        ids[hashval] = ptr;

        /* read the next three lines and ignore them: sequence, header, and quality */
        aline = fgets(line, MAXLINELEN, lfp);
        aline = fgets(line, MAXLINELEN, lfp);
        aline = fgets(line, MAXLINELEN, lfp);
        //fprintf(stderr, "%s", aline);

        nextposition = ftell(lfp);
    }


    /*
     * Now just print all the id lines and their positions
     */

    if (DEBUG) {
        fprintf(stderr, "Bucket sizes\n");

        for (int i = 0; i <= HASHSIZE; i++) {
            struct filehash *ptr;
            ptr = ids[i];
            int counter=0;
            while (ptr != NULL) {
                // fprintf(stdout, "ID: %s Position %ld\n", ptr->id->id, ptr->id->pos);
                counter++;
                ptr = ptr->next;
            }
            fprintf(stderr, "%d\t%d\n", i, counter);
        }
    }

   /* now we want to open output files for left_paired, right_paired, and right_single */

    FILE * left_paired;
    FILE * left_single;
    FILE * right_paired;
    FILE * right_single;

    char *lpfn = catstr(dupstr(left_fn), ".paired.fq");
    char *rpfn = catstr(dupstr(right_fn), ".paired.fq");
    char *lsfn = catstr(dupstr(left_fn), ".single.fq");
    char *rsfn = catstr(dupstr(right_fn), ".single.fq");

    printf("Writing the paired reads to %s and %s.\nWriting the single reads to %s and %s\n", lpfn, rpfn, lsfn, rsfn);

    if ((left_paired = fopen(lpfn, "w")) == NULL ) {
        fprintf(stderr, "Can't open file %s\n", lpfn);
        exit(1);
    }

    if ((left_single = fopen(lsfn, "w")) == NULL) {
        fprintf(stderr, "Can't open file %s\n", lsfn);
        exit(1);
    }
    if ((right_paired = fopen(rpfn, "w")) == NULL) {
        fprintf(stderr, "Can't open file %s\n", rpfn);
        exit(1);
    }

    if ((right_single = fopen(rsfn, "w")) == NULL) {
        fprintf(stderr, "Can't open file %s\n", rsfn);
        exit(1);
    }



    /*
    * Now read the second file, and print out things in common
    */


    FILE *rfp;

    if ((rfp = fopen(right_fn, "r")) == NULL) {
        fprintf(stderr, "Can't open file %s\n", left_fn);
        exit(1);
    }

    nextposition = 0;

    while ((aline = fgets(line, MAXLINELEN, rfp)) != NULL) {

        // make a copy of the current line so we can print it out later.
        char *headerline = dupstr(line);

        line[strcspn(line, " ")] = '\0';

        /* remove the last character, as we did above */

        char lastchar = line[strlen(line)-1];
        if ('1' == lastchar || '2' == lastchar || 'f' == lastchar ||  'r' == lastchar)
            line[strlen(line)-1] = '\0';

        // now see if we have the mate pair
        unsigned hashval = hash(line);
        struct filehash *ptr = ids[hashval];
        long int posn = -1; // -1 is not a valid file position
        while (ptr != NULL) {
            if (strcmp(ptr->id->id, line) == 0) {
                posn = ptr->id->pos;
                ptr->id->printed = true;
            }
            ptr = ptr->next;
        }

        if (posn != -1) {
            // we have a match.
            // lets process the left file
            fseek(lfp, posn, SEEK_SET);
            for (int i=0; i<=3; i++) {
                aline = fgets(line, MAXLINELEN, lfp);
                fprintf(left_paired, "%s", line);
            }
            // now process the right file
            fprintf(right_paired, "%s", headerline);
            for (int i=0; i<=2; i++) {
                aline = fgets(line, MAXLINELEN, rfp);
                fprintf(right_paired, "%s", line);
            }
        }
        else {
            fprintf(right_single, "%s", headerline);
            for (int i=0; i<=2; i++) {
                aline = fgets(line, MAXLINELEN, rfp);
                fprintf(right_single, "%s", line);
            }
        }
    }


    /* all that remains is to print the unprinted singles from the left file */

    for (int i = 0; i <= HASHSIZE; i++) {
        struct filehash *ptr;
        ptr = ids[i];
        while (ptr != NULL) {
            if (! ptr->id->printed) {
                fseek(lfp, ptr->id->pos, SEEK_SET);
                for (int n=0; n<=3; n++) {
                    aline = fgets(line, MAXLINELEN, lfp);
                    fprintf(left_single, "%s", line);
                }
            }
            ptr = ptr->next;
        }
    }


    fclose(lfp);
    fclose(rfp);



    return 0;
}


unsigned hash (char *s) {
    unsigned hashval;

    for (hashval=0; *s != '\0'; s++)
        hashval = *s + 31 * hashval;
    return hashval % HASHSIZE;
}
