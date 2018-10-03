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

uint64_t encode(char *original) {

    size_t length = strlen(original);

    assert(length * 2 == sizeof(uint64_t) * BITS_PER_BYTE);

    uint64_t result = 0;

    for (size_t i = 0; i < length; i++) {
        result = (result << 2) | ((original[i] >> 1) & TWO_BIT_MASK);
    }

    return result;
}

void decode(uint64_t encoded, char *decoded, bool rna_flag) {

    int i = sizeof(uint64_t) * BITS_PER_BYTE / 2;

    for (decoded[i--] = '\0'; i >= 0; i--, encoded >>= 2) {

        unsigned char byte = encoded & TWO_BIT_MASK;

        if (byte == 2) {
            byte = (rna_flag) ? 'U' : 'T';
        } else {
            byte = 'A' | (byte << 1);
        }

        decoded[i] = byte;
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