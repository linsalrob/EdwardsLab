#include <stdio.h>
#include <string.h>
#include "editDistanceLib.h"
//Levenshtein edit distance
//Test
int main() {
	const char *A = "GCGGCAGGAAGGCGCACCCCCCC";
	const char *B = "ACGCGAAGGGAGCCCACCAAAAA";
	printf("%d\n", editDistance(A, B, strlen(A), strlen(B)));
	return 0;
}