#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "FNACharactersLib.h"
#include "queryKmersLib.h"
#include "levenshteinDistanceLib.h"
#include "numberOfLinesLib.h"
//Size of kmers(3 for 3mer, 5 for 5mer, etc. Max 66)
const int KMERSIZE = 5;
const int LINESIZE= 70;
const char *FNANAME1 = "5mers.fna";
//const char *FNANAME2 = "5mers2.fna";
int main() {
	int finalANI = 0;
  int numLines1 = numberOfLines(LINESIZE, FNANAME1) - 1;
  //int numLines2 = numberOfLines(LINESIZE, FNANAME2) - 1;
  char **fna1Characters = fnaCharactersOf(FNANAME1, LINESIZE, numLines1);
  //char **fna2Characters = fnaCharactersOf(FNANAME2, LINESIZE, numLines2);
  int numberOfKmers1 = numLines1*LINESIZE;
  //int numberOfKmers2 = numLines2*LINESIZE;
  char **kmersArrayMalloc1 = queryKmersIntoArray(fna1Characters[0], KMERSIZE, numLines1, LINESIZE);
  free(fna1Characters[0]);
  free(fna1Characters);
  //char **kmersArrayMalloc2 = queryKmersIntoArray(fna2Characters[0], KMERSIZE, numLines2, LINESIZE);
  //free(fna2Characters[0]);
  //free(fna2Characters);
  /*
  int minimumNumberOfKmers;
  if(numberOfKmers1 > numberOfKmers2) {
    minimumNumberOfKmers = numberOfKmers2;
  } else {
    minimumNumberOfKmers = numberOfKmers1;
  }
  for (int k = 0; k < minimumNumberOfKmers; k++) {
		int kmerANI = editDistance(kmersArrayMalloc1[k], kmersArrayMalloc2[k], KMERSIZE, KMERSIZE);
   		finalANI = finalANI + kmerANI;
   	}
    printf("ANI: %d\n", finalANI);
    free(kmersArrayMalloc1);
    free(kmersArrayMalloc2);
    */
    return 0;
}
