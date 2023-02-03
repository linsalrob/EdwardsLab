// convert a fastq file to fasta
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <getopt.h>
using namespace std;

void helpme() {
  fprintf(stderr,
      "fastq2fasta [options] <fastq file> <fasta file>\n"
      "Options:\n"
      "\t-r replace spaces in the fasta header line with '_'\n"
      "\t-n <1/2> add the read number to the header line at the first space\n"
   );
}

int main (int argc, char* argv[]) {

  if ( argc < 3) {
    helpme();
    return 1;
  }

  int replace = 0;
  string n;

  for (;;) {
    switch(getopt(argc, argv, "rn:")) {
      default:
        helpme();
        return 1;
      case -1:
        break;
      case 'r':
        replace = 1;
        continue;
      case 'n':
        n = optarg;
        continue;
    }
    break;
  }
  
  if (optind +2 != argc) {
    helpme();
    return 1;
  }

  char* fqf = argv[optind];
  char* faf = argv[optind+1]; 


  ifstream fastq;
  // streambuf* orig_cin = 0;
  if (strcmp(fqf, "-") != 0) {
	  cout << "reading from " << fqf << '\n';
	  fastq.open(fqf);
	  if (!fastq) return 1;
	  // orig_cin = cin.rdbuf(fastq.rdbuf());
	  cin.rdbuf(fastq.rdbuf());
	  cin.tie(0); // tied to cout by default
  }

  string line;
  ofstream fasta(faf);
  if (!fasta) return 2;
  int c=0;

  while (getline(cin, line))
  {
    //cout << c << " : " << line << '\n';
    if ( c==0 ) {
      line.replace(0, 1, ">");
      if (replace) {
	       for (long unsigned int i = 0; i <= line.length(); i++) {
		       if (line[i] == ' ') {
			       line[i] = '_';
		       }
	       }
      }
      else if (n.size()) {
	      int changed = 0;
	      for (long unsigned int i = 0; i <= line.length(); i++) {
		      if (line[i] == ' ') {
			      line.insert(i++, string("/") + n);
			      changed = 1;
			      break;
		      }
	      }
	      if (changed == 0) {
		      line += string("/") + n;
	      }
      }

      fasta << line << '\n';
    }
    if ( c==1 ) {
      fasta << line << '\n';
    }
    c++;
    if ( c == 4) { c=0 ;}
  }

  if ( c != 0 ) {
    cerr << "ERROR: There appears to be the wrong number of lines in your file!" << endl;
  return 1;
  }

  return 0;
}
