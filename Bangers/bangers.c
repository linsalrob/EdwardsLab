//
// Created by Rob on 17/05/2021.
// Beyond A Nucleotide Generated Evaluation of Relationships between Sequences
// (BANGERS)
// To be used with mash and peas.
//
// Its copyright Rob Edwards and Michael Roach, but good luck if you actually use it.


#include <zlib.h>
#include <stdio.h>
#include <string.h>

#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

//
// We encode the amino acids as a reduced DNA alphabet.
// Note that even though we don't use some of the amino
// acids (e.g. B, O, J, Z, etc, we still include them in
// the array and then we don't need to worry about the 
// conversion
//
 // "A": "GCT", "B": "NNN", "C": "TGT", "D": "AAC",
 // "E": "CAG", "F": "TAT", "G": "GGT", "H": "TAT",
 // "I": "ATT", "J": "NNN", "K": "AAG", "L": "ATT",
 // "M": "ATT", "N": "AAC", "O": "NNN", "P": "CCT",
 // "Q": "CAG", "R": "AAG", "S": "ACT" "T": "ACT",
 // "U": "NNN", "V": "ATT", "W": "TGG", "X": "NNN",
 // "Y": "TAT", "Z": "NNN",
 


int main(int argc, char *argv[])
{

    // declare our array of strings and populate with default values.
    char codons[26][4] = {
	 "GCT", "NNN", "TGT", "AAC", "CAG", "TAT", "GGT",
	 "TAT", "ATT", "NNN", "AAG", "ATT", "ATT", "AAC",
	 "NNN", "CCT", "CAG", "AAG", "ACT", "ACT", "NNN",
	 "ATT", "TGG", "NNN", "TAT", "NNN"
    };

    // set this to 1 to see all the proteins encoded
    int debug = 0;

    // remember gzip reads uncompressed files too!
    gzFile fp;
    kseq_t *seq;
    int l;
    if (argc == 1) {
        fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
        return 1;
    }
    fp = gzopen(argv[1], "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
        printf(">%s\n", seq->name.s);
        // printf("seq: %s\n", seq->seq.s);
	// instantiate an empty string for the encoded version
	char encoded[(strlen(seq->seq.s) * 3) + 1];
	strncpy(encoded, "", sizeof(encoded));
	// enocde the whole protein sequence
	for (int i=0; i < strlen(seq->seq.s); i++) {
		int c = (int) seq->seq.s[i];
		if (c > 64 && c < 91) // upper case amino acids
			c -= 65;
		else if (c > 96 && c < 123) // lower case amino acids
			c -= 97;
		else {
			// something not [A-Za-z]
			printf("Can't parse amino acid: %c\n", seq->seq.s[i]);
			continue;
		}
		if (debug)
			printf("%c\t%d\t%d\t%s\n", seq->seq.s[i], (int) seq->seq.s[i], c, codons[c]);
		
		strcat(encoded, codons[c]);

	}
	printf("%s\n", encoded);

    }
    kseq_destroy(seq); 
    gzclose(fp);
    return 0;
}
