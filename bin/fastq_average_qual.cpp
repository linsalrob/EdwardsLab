// convert a fastq file to fasta
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
using namespace std;

int main (int argc, char* argv[]) {

  if ( argc < 2) {
    cerr << "Usage: " << argv[0] << " <fastq file (use - to read from STDIN)>\n";
    return 1;
  }


  ifstream fastq;
  streambuf* orig_cin = 0;
  if (strcmp(argv[1], "-") != 0) {
    cout << "reading from " << argv[1] << '\n';
    fastq.open(argv[1]);
    if (!fastq) return 1;
    orig_cin = cin.rdbuf(fastq.rdbuf());
    cin.tie(0); // tied to cout by default
  }

  string line;
  int c=0;
  long total=0;
  long n=0;
  while (getline(cin, line))
  {
    //cout << c << " : " << line << '\n';
    if ( c==3 ) {
	for (char const &b: line) {
		total += (int) b;
		n++;
	}
    }

    c++;
    if ( c == 4) { c=0 ;}
  }

  if ( c != 0 ) {
	  cerr << "ERROR: There appears to be the wrong number of lines in your file!" << endl;
	return 1;
  }

  cout << "Bases: " << n << " Total quality: " << total << " Average quality: " << total/n << endl;

  return 0;
}
