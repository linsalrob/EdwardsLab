#include <stdio.h>
#include <string.h>
//Levenshtein edit distance
int editDistance(const char* S, const char* T, int sLength, int tLength) {
	int minimum, first, second, third, conditional;
	//Maximum
	if(sLength == 0) {
		return tLength;
	}
	if(tLength == 0) {
		return sLength;
	}
	if(S[sLength-1] != T[tLength - 1]) {
		conditional = 1;
	} else {
		conditional = 0;
	}

	//Minimum
	first = editDistance(S, T, sLength - 1, tLength) + 1;
	second = editDistance(S, T, sLength, tLength - 1) + 1;
	third = editDistance(S, T, sLength - 1, tLength - 1) + conditional;

	if(first > second) {
		minimum = second;
	} else if (second > third) {
		minimum = third;
	} else {
		minimum = first;
	}
	return minimum;
}

//test
int main() {
	const char *A = "GCGGCAGGAAGGCGCACCCCCCC";
	const char *B = "ACGCGAAGGGAGCCCACCAAAAA";

	printf("%d\n", editDistance(A, B, strlen(A), strlen(B)));

	return 0;
}