#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "numberOfLinesLib.h"

int numberOfLines(const int LINESIZE, const char *FNANAME) {
    FILE *fnaPointer;
    char singleCharLine[256];
    fnaPointer = fopen(FNANAME, "r");
    if(fnaPointer == NULL) {
        printf("%s", "Error: fnaPointer is null\n");
        exit(EXIT_FAILURE);
    }
    int numberOfLines = 0;
    while (fgets(singleCharLine, 256, fnaPointer) != NULL) {
        numberOfLines++;
    }
    if(numberOfLines == 0) {
        printf("%s", "Error: No Lines in .fna file\n");
        exit(EXIT_FAILURE);
    }
    fclose(fnaPointer);
    return numberOfLines;
}