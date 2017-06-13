#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "queryKmersLib3.h"

char* queryKmersIntoArray(const char* dnaChars, int kmerNum) {
	char** kmersArrayMalloc = malloc(strlen(dnaChars)*sizeof(char*));
	if(kmersArrayMalloc == NULL) {
		printf("%s", "kmersArrayMallocFailed");
	}
	//char kmersArrayH[50000];
	char buff[100];
	for(int i = 0; i < strlen(dnaChars) - kmerNum + 1; i++) {
		//printf("%s",kmersArrayMalloc);
		//memcpy(buff, dnaChars[i], kmerNum); can do normal to pointer
		//printf("%c",buff);
		//kmersArrayMalloc[i] = strdup(dnaChars + i, kmerNum);
		//kmersArrayMalloc[i] = *strndup(dnaChars + i, kmerNum);
		//char* x = strndup(dnaChars + i, kmerNum);
		//strcpy(kmersArrayMalloc[i], x);
		//char sbuff[100];
		//kmersArrayMalloc[i] = 
		//memcpy(kmersArrayMalloc[i], dnaChars + i, kmerNum);
		//printf("%c",sbuff);
		memcpy(buff, &dnaChars[i], kmerNum);
		//*(kmersArrayMalloc + i) = buff;
		kmersArrayMalloc[i] = buff;
		//*(kmersArrayMalloc + i) = buff;
		//strcat(kmersArrayMalloc[i], buff);
		//printf("%s", buff);

		//strncpy(kmersArrayMalloc, &dnaChars[i], kmerNum);
	}
	printf("%s", *kmersArrayMalloc);
		//printf("%s",kmersArrayMalloc);
	//char* kmersArrayHeap = strdup(kmersArrayMalloc);
	//printf("%s", kmersArrayHeap);
	return (char*) kmersArrayMalloc;
}

