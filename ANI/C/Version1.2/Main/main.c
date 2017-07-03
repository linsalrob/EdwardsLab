#include <stdio.h>
#include <dirent.h>
#include <string.h>
#include <stdlib.h>

#include "numberOfLinesLib.h"
#include "FNACharactersLib.h"
#include "levenshteinDistanceLib.h"

const int DATABUFFER2[4096];
const int LINESIZE = 70;
const int KMERSIZE = 15;
const int fnaFilesDirectoryStringLength = 41 + 2;
const char fnaFilesDirectory[fnaFilesDirectoryStringLength] = {"."};
const int minimumFnaFileStringLength = 4;
const int fnaNameBufferSize = 4;
const int amountOfFNAFiles = 6;
const int nameLength = 200;

int main() {
    int fnaFileCount = 0;
    DIR *fnaFilesDirectoryPointer;
    struct dirent *direntFnaDirectory;
    fnaFilesDirectoryPointer = opendir(fnaFilesDirectory);
    char checkIfFNA[fnaNameBufferSize];
    char fnaString[fnaNameBufferSize] = "fna";
    int fnaFilesDirectoryStringLengthMinusOne = minimumFnaFileStringLength - 1;
    char *x[amountOfFNAFiles];
    char** fnaArrayMalloc = (char**) malloc(amountOfFNAFiles*sizeof(char*));
    if(fnaArrayMalloc == NULL) { 
        printf("%s", "Error: kmersArrayMalloc\n");
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < amountOfFNAFiles; i++) {
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
                    memcpy(fnaArrayMalloc[fnaFileCount], &direntFnaDirectory->d_name, strlenString);
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
      for(int i=0; i < amountOfFNAFiles; i++) {
        for(int j=i; j < amountOfFNAFiles; j++) {
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
                } 
                if(booleanVal1 == 0) {
                    int kmerANI = editDistance(tempArr, tempArr2, KMERSIZE, KMERSIZE);
                    finalANI = finalANI + kmerANI;
                    kmerCount++;
                } 
            }
            printf("\n %s %s ANI: %d\n", fnaArrayMalloc[i], fnaArrayMalloc[j], finalANI);
            finalANI = 0;
            }
        }
}