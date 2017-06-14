#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "queryKmersLib.h"

char** queryKmersIntoArray(const char* dnaChars, int kmerSize, int numKmers) {

	char** kmersArrayMalloc = (char**) malloc(numKmers*sizeof(char*));
	if(kmersArrayMalloc == NULL) { printf("%s", "1");}
	for(int i = 0; i < numKmers; i++) {
      kmersArrayMalloc[i] = (char *) malloc(sizeof(char)*kmerSize);
   }

	for(int i = 0; i < strlen(dnaChars) - kmerSize + 1; i++) {
		printf("%s", &dnaChars[i]);
		printf("%s", " ");
		memcpy(kmersArrayMalloc[i], &dnaChars[i], kmerSize);
	}
	return (char**) kmersArrayMalloc;
}