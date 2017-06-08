#include <stdio.h>
#include <string.h>
#include "editDistanceLib.h"
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
	//Conditional 1 or 0
	if(S[sLength-1] != T[tLength - 1]) {
		conditional = 1;
	} else {
		conditional = 0;
	}
	//Minimum
	first = editDistance(S, T, sLength - 1, tLength) + 1;
	second = editDistance(S, T, sLength, tLength - 1) + 1;
	third = editDistance(S, T, sLength - 1, tLength - 1) + conditional;
	minimum = first;
	if(first > second) {
		minimum = second;
	} 
	if (second > third) {
		minimum = third;
	}
	return minimum;
}