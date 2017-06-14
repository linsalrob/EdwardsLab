#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "FNACharactersLib.h"

const int NUMLINES = 70000;
const int LINESIZE= 70;
const int BUFF = 100;

char** fnaCharactersOf(const char* fnaFileNameAndLocation) {
    FILE *fnaPointer;
    const char* fnaName = fnaFileNameAndLocation;
    char singleCharLine[BUFF], firstLineOfChars[BUFF];
    fnaPointer = fopen(fnaName, "r");
    fgets(firstLineOfChars, BUFF, fnaPointer);
    
    char **fnaArrayMalloc = (char**) malloc(NUMLINES*sizeof(char*));
    if(fnaArrayMalloc == NULL) { printf("%s", "1");}

    for(int i = 0; i < NUMLINES; i++) {
        fnaArrayMalloc[i] = (char *) malloc(sizeof(char)*LINESIZE);
   }
   
   int j = 0;
    while (fgets(singleCharLine, BUFF, fnaPointer) != NULL) {
        memcpy(fnaArrayMalloc[j], &singleCharLine, LINESIZE);
        j++;
    }

    char **fnaCharMalloc = (char**) malloc(sizeof(char));
    if(fnaCharMalloc == NULL) { printf("%s", "1");}
    fnaArrayMalloc[0] = (char *) malloc(sizeof(char)*LINESIZE*NUMLINES);


    //for (int k = 0; k < NUMLINES; k++) {
        //printf("%s", fnaArrayMalloc[k]);
        //memcpy(fnaCharMalloc[0], &fnaArrayMalloc[k], LINESIZE);
        //memcpy(fnaCharMalloc[0] + LINESIZE, &fnaArrayMalloc[k], LINESIZE);
        //memcpy(fnaCharMalloc[LINESIZE] + LINESIZE, &fnaArrayMalloc[k], LINESIZE);
        //strcpy(fnaCharMalloc[0], fnaArrayMalloc[0]);
    //}

    printf("%s",fnaCharMalloc[0]);

    fclose(fnaPointer);
    return (char**) fnaArrayMalloc;
}