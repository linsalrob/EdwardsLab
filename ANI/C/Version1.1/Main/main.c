#include <stdio.h>
#include <string.h>
#import "numberOfLinesLib.h"
#import "FNACharactersLib.h"
#import "levenshteinDistanceLib.h"
const int DATABUFFER2[4096];
const int LINESIZE = 70;
const int KMERSIZE = 100;
const char *FNANAME1 = "NC_021215.fna";
const char *FNANAME2 = "NC_022886.fna";

int main() {
    const int kmerSizeMod = KMERSIZE+1;
    int finalANI = 0;
    char *ptr, *ptr2;
    char tempArr[kmerSizeMod], tempArr2[kmerSizeMod];
    FILE *fnaPointer;
    int numLines1 = numberOfLines(LINESIZE, FNANAME1);
    int numLines2 = numberOfLines(LINESIZE, FNANAME2);
    char **fna1Characters = fnaCharactersOf(FNANAME1, LINESIZE, numLines1); 
    char **fna2Characters = fnaCharactersOf(FNANAME2, LINESIZE, numLines2);
    int minimumLines;
    ptr=fna1Characters[0];
    ptr2=fna2Characters[0];
    int chr = 0;
    int chr2 = 0;
    int kmerCount = 0;
    if(numLines1 < numLines2) {
        minimumLines = numLines1;
    } else {
        minimumLines = numLines2;
    }
    const int numChars = minimumLines*LINESIZE;
    const int numKmers = (numLines1-1)*LINESIZE-KMERSIZE;
    for(int i=0; i<numChars; i++) {
        for(int j=0; j<kmerSizeMod; j++) {
            if(j == 0) {
                chr = i;
                chr2 = i;
            }
            if(j==kmerSizeMod-1) {
                tempArr[j] = '\0';
            } else {
                tempArr[j] = ptr[chr++];
            }
            if(j==kmerSizeMod-1) {
                tempArr2[j] = '\0';
            } else {
                tempArr2[j] = ptr2[chr2++];
            }
        } 
        if(strlen(tempArr) != KMERSIZE && strlen(tempArr2) != KMERSIZE) {break;}
        int kmerANI = editDistance(tempArr, tempArr2, KMERSIZE, KMERSIZE);
        finalANI = finalANI + kmerANI;
        //printf("\nTempArr1: %s\n", tempArr);
        //printf("TempArr2: %s\n", tempArr2);
        kmerCount++;
        //printf("%d\n", kmerANI);
    }
    printf("\n ANI: %d\n", finalANI);
}