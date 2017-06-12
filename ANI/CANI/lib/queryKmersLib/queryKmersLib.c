#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "queryKmersLib.h"

char* queryKmersIntoArray(const char* dnaChars, int kmerNum) {
	char *kmersArray[strlen(dnaChars)]; 
	for(int i = 0; i < strlen(dnaChars) - kmerNum + 1; i++) {
		kmersArray[i] = strndup(dnaChars + i, kmerNum);
	}
	printf("%s\n", kmersArray[i]);
	return *kmersArray;
}

