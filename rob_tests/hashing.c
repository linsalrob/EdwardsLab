/*
 * Calculate some hashes for sequences, and also try and generate O(n) code
 * to calculate hash for next hash in a string
 *
 * This should be O(n) complexity where n = length of sequence
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

// alpha for our hash function
#define alpha 5

long double hash(char* seq) {
	// calculate the hash function for a single string
	
	long double h = 0;
	for (int i=0; i<strlen(seq); i++)
		h += (pow(alpha, ((int)strlen(seq) - i))) * seq[i];
	
	return h;
}

long double next_hash(char* seq, int p, int k, long double last_hash) {
	
	return (long double) alpha * (last_hash - (pow((long double) alpha, (long double) k) * (long double) seq[p-1]) + (long double) seq[p+k-1]);
}

unsigned long long ahash(char* seq) {
	unsigned long long h = 0;
	int c = 0;
	while ((c=*seq++)) {
		unsigned long long nh = alpha * (h + c);
		if ((LLONG_MAX - h) < nh)
			fprintf(stderr, "\033[91mWARNING: %lld will be exceeded!!\033[0m\n", LLONG_MAX);
		h = alpha * (h + c);
	}
	return h;
}



int main(int argc, char *argv[]) {
	char* s = "abc";
	char* t = "bcd";
	char* u = "abcd";

	long double hs = hash(s);
	long double ht = hash(t);
	long double hu = next_hash(u, 1, (int) strlen(s), hs);
	long double diff = hu - ht;
	printf("Hash of %s is %Lf. Hash of %s is %Lf. Hash of substring is %Lf. Delta is %Lf\n", s, hs, t, ht, hu, diff);

	s = "abcdefghijkl";
	t = "bcdefghijklm";
	u = "abcdefghijklm";

	hs = hash(s);
	ht = hash(t);
	hu = next_hash(u, 1, (int) strlen(s), hs);
	diff = hu - ht;
	printf("Hash of %s is %Lf. Hash of %s is %Lf. Hash of substring is %Lf. Delta is %Lf\n", s, hs, t, ht, hu, diff);

	unsigned long long hhs = ahash(s);
	unsigned long long hht = ahash(t);
	printf("aHash of %s is %lld. aHash of %s id %lld\n", s, hhs, t, hht);

}
	


