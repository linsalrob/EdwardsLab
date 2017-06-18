#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "FNACharactersLib.h"

const int DATABUFFER = 100;
char** fnaCharactersOf(const char* fnaFileNameAndLocation, const int LINESIZE, const int numberOfLines) {
    FILE *fnaPointer;
    const char* fnaName = fnaFileNameAndLocation;
    char singleCharLine[DATABUFFER], firstLineOfChars[DATABUFFER];
    fnaPointer = fopen(fnaName, "r");
    if(fnaPointer == NULL) {
        printf("%s", "Error: fnaPointer is null\n");
        exit(EXIT_FAILURE);
    }//Get first line of .fna file(doesn't include respective ACGT)
    fgets(firstLineOfChars, DATABUFFER, fnaPointer);
    char **fnaArrayMalloc = (char**) malloc(numberOfLines*sizeof(char*));
    if(fnaArrayMalloc == NULL) {
        printf("%s", "Error: fnaArrayMalloc\n");
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < numberOfLines; i++) {
        fnaArrayMalloc[i] = (char *) malloc(sizeof(char)*LINESIZE);
   }
   
   int j = 0;
    while (fgets(singleCharLine, DATABUFFER, fnaPointer) != NULL) {
        memcpy(fnaArrayMalloc[j], &singleCharLine, LINESIZE);
        j++;
    }
    fclose(fnaPointer);
    char **fnaCharMalloc = (char**) malloc(sizeof(char));
    if(fnaCharMalloc == NULL) {
        printf("%s", "Error: fnaCharMalloc\n");
        exit(EXIT_FAILURE);
    }
    fnaCharMalloc[0] = (char *) malloc(sizeof(char)*LINESIZE*numberOfLines);
    for (int k = 0; k < numberOfLines; k++) {
        memcpy(fnaCharMalloc[0] + LINESIZE*k, fnaArrayMalloc[k], LINESIZE);
    }
    for(int l = 0; l < numberOfLines; l++) {
        free(fnaArrayMalloc[l]);
    }
    free(fnaArrayMalloc);
    
    return (char**) fnaCharMalloc;
}