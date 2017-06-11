#include <stdio.h>
#include <string.h>
#include "FNACharactersLib.h"

char* fnaCharactersOf(const char* fnaFileNameAndLocation) {
    FILE *fnaPointer;
    const char* fnaName = fnaFileNameAndLocation;
    char singleCharLine[100], firstLineOfChars[100], dnaChars[8000000];
    fnaPointer = fopen(fnaName, "r");
    //if (fnaPointer == NULL){ return 1; }
    fgets(firstLineOfChars, 100, fnaPointer);
    while (fgets(singleCharLine, 100, fnaPointer) != NULL) { 
    	strcat(dnaChars, singleCharLine);
    }
    fclose(fnaPointer);
    char* dnaCharsHeap = strdup(dnaChars);//or memcopy(sizeOf)
    return dnaCharsHeap;
    return 0;
}