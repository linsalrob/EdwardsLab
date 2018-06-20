// convert a fastq file to fasta
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
using namespace std;

int main (int argc, char* argv[]) {

  if ( argc < 3) {
    cerr << "Usage: " << argv[0] << " <fastq file (use - to read from STDIN)> <fasta file>\n";
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
  ofstream fasta(argv[2]);
  if (!fasta) return 2;
  int c=0;
  
  while (getline(cin, line))
  {
    //cout << c << " : " << line << '\n';
    if ( c==0 ) {
      line.replace(0, 1, ">");
      fasta << line << '\n';
    }
    if ( c==1 ) {
      fasta << line << '\n';
    }
    c++;
    if ( c == 4) { c=0 ;}
  }

  return 0;
}
