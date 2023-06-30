/*
 * Find primers whereever they are in the sequence. We hash the primers and then look through the sequence to see if we have that hash
 *
 * This should be O(n) complexity where n = length of sequence
 */


#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "kseq.h"
#include "version.h"
#include <float.h>

// alpha for our hash function
#define alpha 33

long double hash(char* seq) {
	// calculate the hash function for a single string
	
	long double h = 0;
	for (int i=0; i<strlen(seq); i++) {
		long double d = (pow(alpha, ((int)strlen(seq) - i))) * seq[i];
		double sls = ((int)strlen(seq) - i);
		double psls = pow(alpha, sls);
		h += d;
		double dd = psls * seq[i];
		if (dd != d)
			fprintf(stderr, "\033[91mWARNING: %Lf is not %f\033[0m\n", d, dd);

		fprintf(stderr, "%s : %c : (%d^%f) = %f. (%f * %d) = %Lf. Hash: %Lf\n ", seq, seq[i], alpha, sls, psls, psls, seq[i], d, h);   
	}
	fprintf(stderr, "%s : %Lf\n", seq, h); 
	return h;
}

long double next_hash(char* seq, int p, int k, long double last_hash) {
	// For a string of length k starting at p 
	// given the previous string starting at p-1 has hash last_hash
	// calculate the next hash
	//
	
	double d = pow(alpha, k);
	
	fprintf(stderr, "%d * (%Lf - (%f * %d) + %d)\n", alpha, last_hash, d, seq[p-1], seq[p+k-1]);  

	long double h = (long double) alpha * (last_hash - (pow((long double) alpha, (long double) k) * (long double) seq[p-1]) + (long double) seq[p+k-1]);
	fprintf(stderr, "hash at %d : %Lf\n", p, h);
	return h;
}

unsigned long long int uli_hash(char* seq) {
	
	unsigned long long int h = 0;
	// for (int i=0; i<strlen(seq); i++) 
	for (int i=strlen(seq)-1; i>=0; i--) 
		h += (unsigned long long int) (pow(alpha, ((int)strlen(seq) - i))) * seq[i];
	
	fprintf(stderr, "%s : %llu\n", seq, h); 
	return h;
}
unsigned long long int another_hash(char* seq, int p, int k, unsigned long long int last_hash) {
	unsigned long long int h = alpha * (last_hash - ((unsigned long long int) pow(alpha, k) * seq[p-1]) + seq[p+k-1]);
	fprintf(stderr, "ULI hash at %d : %llu\n", p, h);
	return h;
}


unsigned long dbj_hash(char *seq)
{
	unsigned long hash = 0;
	int c = 0;

	while ((c = *seq++)) 
		hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

	fprintf(stderr, "DBJ hash : %lu\n", hash);
	return hash;
}


int main(int argc, char *argv[]) {
	/*
	   char* s = "abcdefghij";
	   char* t = "bcdefghijk";
	   char* u = "abcdefghijk";
	char* s = "abcdefghijklmn";
	char* t = "bcdefghijklmno";
	char* u = "abcdefghijklmno";
	*/

	char* s = "abc";
	char* t = "bcd";
	char* u = "abcd";

	/*
	long double hs = hash(s);
	long double ht = hash(t);
	long double hu = next_hash(u, 1, (int) strlen(s), hs);
	long double diff = hu - ht;
	fprintf(stderr, "Delta is %Lf\n", diff);
	*/
	unsigned long long int hs = uli_hash(s);
	unsigned long dbjh = dbj_hash(s);
	unsigned long long int ht = uli_hash(t);
	unsigned long long int hu = another_hash(u, 1, (int) strlen(s), hs);
	unsigned long long int diff = ht - hu;
	fprintf(stderr, "Delta is %llu\n", diff);

}
	


