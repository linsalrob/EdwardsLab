#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "queryKmersLib.h"
// 2 <= kmerSize <= 66
//67+ will return kmers of 1000+ characters for some reason(lineSize = 70?)
char** queryKmersIntoArray(const char* dnaChars, const int KMERSIZE, const int numberOfLines, const int LINESIZE) {
	const int numberOfKmers = numberOfLines*LINESIZE;
	char** kmersArrayMalloc = (char**) malloc(numberOfKmers*sizeof(char*));
	if(kmersArrayMalloc == NULL) { 
        printf("%s", "Error: kmersArrayMalloc\n");
        exit(EXIT_FAILURE);
	}
	for(int i = 0; i < numberOfKmers; i++) {
      kmersArrayMalloc[i] = (char *) malloc(sizeof(char)*KMERSIZE*10);
      	if(kmersArrayMalloc[i] == NULL) { 
        printf("%s", "Error: kmersArrayMalloc\n");
        exit(EXIT_FAILURE);
	}
   }
	for(int i = 0; i < strlen(dnaChars) - KMERSIZE + 1; i++) {
		memcpy(kmersArrayMalloc[i], &dnaChars[i], KMERSIZE);
	}
	return (char**) kmersArrayMalloc;
	
}