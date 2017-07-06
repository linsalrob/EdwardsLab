#include <stdio.h>
#include <dirent.h>
#include <string.h>
#include <stdlib.h>

#include "numberOfLinesLib.h"
#include "FNACharactersLib.h"
#include "levenshteinDistanceLib.h"

const int LINESIZE = 256;
const int KMERSIZE = 2;
const int fnaFilesDirectoryStringLengthBuffer = 100 + 2;
const char fnaFilesDirectory[fnaFilesDirectoryStringLengthBuffer] = {"/Users/jon/Desktop/Version133/Main/fna/"};
const int minimumFnaFileStringLength = 4;
const int fnaNameBufferSize = 4;
const int amountOfFNAFilesBuffer = 2200;
const int nameLength = 100;

int main() {
    int fnaFileCount = 0;
    DIR *fnaFilesDirectoryPointer;
    struct dirent *direntFnaDirectory;
    fnaFilesDirectoryPointer = opendir(fnaFilesDirectory);
    char checkIfFNA[fnaNameBufferSize];
    char fnaString[fnaNameBufferSize] = "fna";
    int fnaFilesDirectoryStringLengthMinusOne = minimumFnaFileStringLength - 1;
    //char *x[amountOfFNAFilesBuffer];
    char** fnaArrayMalloc = (char**) malloc(amountOfFNAFilesBuffer*sizeof(char*));
    if(fnaArrayMalloc == NULL) { 
        printf("%s", "Error: kmersArrayMalloc\n");
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < amountOfFNAFilesBuffer; i++) {
        fnaArrayMalloc[i] = (char *) malloc(sizeof(char)*nameLength);
        if(fnaArrayMalloc[i] == NULL) { 
            printf("%s", "Error: kmersArrayMalloc\n");
            exit(EXIT_FAILURE);
        }
   }
    if(fnaFilesDirectoryPointer) {
        while ((direntFnaDirectory = readdir(fnaFilesDirectoryPointer)) != NULL) {
            int strlenString = strlen(direntFnaDirectory->d_name);
            if( strlenString > minimumFnaFileStringLength) {
                memcpy(checkIfFNA, &direntFnaDirectory->d_name[strlenString - fnaFilesDirectoryStringLengthMinusOne], fnaFilesDirectoryStringLengthMinusOne);
                if (strcmp(checkIfFNA,fnaString) == 0) {
                    int strlenString2 = strlen(fnaFilesDirectory);
                    memcpy(fnaArrayMalloc[fnaFileCount], fnaFilesDirectory, strlenString2);
                    memcpy(fnaArrayMalloc[fnaFileCount] + strlenString2, &direntFnaDirectory->d_name, strlenString);
                    fnaFileCount++;
                }
            }
        }
        closedir(fnaFilesDirectoryPointer);
    } else {
        printf("%s", "Error: fnaFilesDirectory[fnaFilesDirectoryStringLength] is null\n");
        exit(EXIT_FAILURE);
    }
    const int kmerSizeMod = KMERSIZE+1;
    int finalANI = 0;
    char *ptr, *ptr2;
    char tempArr[kmerSizeMod], tempArr2[kmerSizeMod];
    FILE *fnaPointer;
      for(int i=0; i < fnaFileCount; i++) {
        for(int j=i; j < fnaFileCount; j++) {
            int numLines1 = numberOfLines(LINESIZE, fnaArrayMalloc[i]) - 1;          
            int numLines2 = numberOfLines(LINESIZE, fnaArrayMalloc[j]) - 1;
            char **fna1Characters = fnaCharactersOf(fnaArrayMalloc[i], LINESIZE, numLines1);
            char **fna2Characters = fnaCharactersOf(fnaArrayMalloc[j], LINESIZE, numLines2);
            int minimumLines;
            int maximumLines;
            int chr = 0;
            int chr2 = 0;
            int kmerCount = 0;
            ptr=fna1Characters[0];
            ptr2=fna2Characters[0];
            if(numLines1 < numLines2) {
                minimumLines = numLines1;
                maximumLines = numLines2;
            } else {
                minimumLines = numLines2;
                maximumLines = numLines1;
            }
            const int numCharsMin = minimumLines*LINESIZE;
            const int numCharsMax = maximumLines*LINESIZE;
            const int numKmers = (numLines1-1)*LINESIZE-KMERSIZE;
            int booleanVal1 = 0;
            
            for(int i=0; i<numCharsMin; i++) {
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
                /*
                if(strlen(tempArr) != KMERSIZE) {
                    if(booleanVal1 == 1) {
                    } else {
                    //printf("\n((numCharsMax-i)/KMERSIZE): %d", ((numCharsMax - i)/KMERSIZE));
                    finalANI = finalANI + (numCharsMax - i)/KMERSIZE;
                    booleanVal1 = 1;
                    }
                }
                if(strlen(tempArr2) != KMERSIZE) {
                    if(booleanVal1 == 1) {
                    } else {
                    //printf("\n((numCharsMax-i)/KMERSIZE): %d", ((numCharsMax - i)/KMERSIZE));
                    finalANI = finalANI + (numCharsMax - i)/KMERSIZE;
                    booleanVal1 = 1;
                    }
                } */
                if(i == 100) {
                    //
                    break;
                }
                if(booleanVal1 == 0) {
                    int kmerANI = editDistance(tempArr, tempArr2, KMERSIZE, KMERSIZE);
                    finalANI = finalANI + kmerANI;
                    kmerCount++;
                    //
                }
            }
                //}
            //}
            double total = numCharsMax/KMERSIZE;
            double sim = (finalANI/total);
            printf("\n %s %s ANI: %d\n", fnaArrayMalloc[i], fnaArrayMalloc[j], finalANI);
            //printf("\n %s %s Percent Simularity: %f\n", fnaArrayMalloc[i], fnaArrayMalloc[j], sim);
            free(fna1Characters);
            free(fna2Characters);
            finalANI = 0;
            }
        }
}