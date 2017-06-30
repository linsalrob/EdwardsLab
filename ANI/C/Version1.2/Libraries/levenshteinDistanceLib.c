#include <stdio.h>
int editDistance(const char* S, const char* T, int kmerLength) {
    int matrix[kmerLength + 1][kmerLength + 1];
    for (int i = 0; i <= kmerLength; i++) {
        matrix[i][0] = i;
        matrix[0][i] = i;
    }
    for (int j = 1; j <= kmerLength; j++) {
        for (int k = 1; k <= kmerLength; k++) {
            if (S[j-1] == T[k-1]) 
                matrix[j][k] = matrix[j-1][k-1];
            else 
                matrix[j][k] = (matrix[j][k-1] + 1) > (matrix[j-1][k] + 1) ? ((matrix[j-1][k-1] + 1) > (matrix[j-1][k] + 1) ? (matrix[j-1][k] + 1) : (matrix[j-1][k-1] + 1)) : ((matrix[j-1][k-1] + 1) > (matrix[j][k-1] + 1) ? (matrix[j][k-1] + 1) : (matrix[j-1][k-1] + 1));
        }
    }
    return matrix[kmerLength][kmerLength];
}