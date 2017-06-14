#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "FNACharactersLib.h"
#include "queryKmersLib.h"
#include "editDistanceLib.h"

//Max 100,000(100000) for number or fault will occur
const int NUMKMERS = 60;//Number of kmers(might change to dynamic) min 100
//Max ? for speed of edit distance
const int KMERSIZE = 5;//Size of kmers(3 for 3mer, 5 for 5mer, etc.)
//
const int NUMLINES2 = 70000;
//
const char *fnaName = "NC_010468.fna";
//
const char *fnaName2 = "5mers.fna";
//
int main() {
	//int finalANI = 0;
	//FNACharactersLib
	//1
  //char* fnaName = "5mers.fna";
  char **fnaArrayMalloc = fnaCharactersOf(fnaName); 
  /*
  for (int i = 0; i < NUMLINES2; i++) {
      printf("Value of names[%d] = %s\n", i, fnaArrayMalloc[i]);
   }
   */
    //2
    //char *fna2 = fnaCharactersOf(fnaName2);
    //queryKmersLib
    //1
   /*
    char **kmersArrayMalloc = queryKmersIntoArray(fna1, KMERSIZE, NUMKMERS);
    for (int i = 0; i < NUMKMERS; i++) {
      printf("Value of names[%d] = %s\n", i, kmersArrayMalloc[i]);
   }
    free(fna1);
    free(kmersArrayMalloc);
    //2
    char **kmersArrayMalloc2 = queryKmersIntoArray(fna2, KMERSIZE, NUMKMERS);
    for (int j = 0; j < NUMKMERS; j++) {
      printf("Value of names[%d] = %s\n", j, kmersArrayMalloc2[j]);
   }
   free(fna2);
   */
   //editDistanceLib
   //1 and 2
   /*
   for (int k = 0; k < NUMKMERS; k++) {
		int kmerANI = editDistance(kmersArrayMalloc[k], kmersArrayMalloc2[k], KMERSIZE, KMERSIZE);
   		finalANI = finalANI + kmerANI;
   	}
    free(kmersArrayMalloc);
    free(kmersArrayMalloc);
    printf("ANI: %d", finalANI);
    */
    return 0;
}