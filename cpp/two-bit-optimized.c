// Two bit DNA encoding
// Note, this solution comes from https://stackoverflow.com/questions/39242932/how-to-encode-char-in-2-bits
// Created by redwards on 10/1/18.
//

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>

#define TWO_BIT_MASK (3)
#define BITS_PER_BYTE (8)
#define BIG_ENOUGH (1024)

char encodeTable [] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 'B', 1, 'D', 'E', 'F', 3, 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 2, 2, 'V', 'W', 'X', 'Y', 'Z'};
char decodeTableDNA [] = {'A', 'C', 'T', 'G'};
char decodeTableRNA [] = {'A', 'C', 'U', 'G'};

uint64_t inline encode(char *original) {
    uint64_t result = ((original[0] >> 1) & TWO_BIT_MASK);
    for (size_t i = 1; i < 32; i++) {
        result = (result << 2) | encodeTable[original[i]];
    }

    return result;
}


void inline decode(uint64_t encoded, char *decoded, bool rna_flag) {

    int i = sizeof(uint64_t) * BITS_PER_BYTE / 2;

    if (rna_flag) {
        for (decoded[i--] = '\0'; i >= 0; i--, encoded >>= 2) {

            unsigned char byte = encoded & TWO_BIT_MASK;

            byte = decodeTableRNA[byte];

            decoded[i] = byte;
        }
    }
    else {
        for (decoded[i--] = '\0'; i >= 0; i--, encoded >>= 2) {

            unsigned char byte = encoded & TWO_BIT_MASK;

            byte = decodeTableDNA[byte];

            decoded[i] = byte;
        }
    }
}

int main() {
    char *segment = "GCCGTGCTAAGCGTAACAACTTCAAATCCGCG";

    printf("%s\n", segment);

    uint64_t binary = encode(segment);

    printf("%llu\n", binary);

    char string[BIG_ENOUGH];

    decode(binary, string, false);

    printf("%s\n", string);

    return 0;
}
