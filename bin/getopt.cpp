#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <getopt.h>
using namespace std;

int main (int argc, char* argv[]) {
	if ( argc < 3 ) {
		cerr << "Please add some options\n";
		return 1;
	}

	int replace = 0;
	for (;;) {
		switch(getopt(argc, argv, "r")) {
			default:
				cerr << " fount " << optarg << "\n";
			case -1:
				cerr << "Found -1\n";
				break;
			case 'r':
				cout << "Flag r set";
				replace = 1;
				continue;
		}
		break;
	}
	cerr << "argc: " << argc << "\n";
	cerr << "Optind: " << optind << "\n";
	if (optind +2 != argc) {
		cerr << " <fastq file> <fasta file>";
		return 1;
	}
	char* fqf = argv[optind];
	char* faf = argv[optind+1]; 
	cout << "Fastq: " << fqf << " Fasta: " << faf << "\n";
}




