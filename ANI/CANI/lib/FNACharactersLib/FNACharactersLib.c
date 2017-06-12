#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "FNACharactersLib.h"

char* fnaCharactersOf(const char* fnaFileNameAndLocation) {
    FILE *fnaPointer;
    const char* fnaName = fnaFileNameAndLocation;
    char singleCharLine[100], firstLineOfChars[100], dnaChars[8000000];
    fnaPointer = fopen(fnaName, "r");
    fgets(firstLineOfChars, 100, fnaPointer);
    char x[100];
    while (fgets(singleCharLine, 100, fnaPointer) != NULL) {
        for(int i = 0; i < strlen(singleCharLine); i++) {
            if(isspace((char)singleCharLine[i]) != 1) {
                strcat(dnaChars, &singleCharLine[i]);
                //printf("%c\n", (char)singleCharLine[i]);
            }
            //printf("%s\n",dnaChars);
        }
    }
    fclose(fnaPointer);
    char* dnaCharsHeap = strdup(dnaChars);//or memcopy(sizeOf)
    return dnaCharsHeap;
}