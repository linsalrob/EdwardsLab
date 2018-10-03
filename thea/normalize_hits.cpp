//
// Created by redwards on 10/3/18.
//
// Normalize the rapsearch hits. This just sums the counts and then prints the numbers
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <map>
using namespace std;

int main (int argc, char* argv[]) {

    if ( argc < 3) {
        cerr << "Usage: " << argv[0] << " <file to normalize> <directory to write to>\n";
        return 1;
    }

    // split the input file name so we can get an output filename
    string ofname = (string) argv[1];
    auto pos = ofname.rfind('/');
    if (pos != string::npos)
        ofname = ofname.substr(pos+1, string::npos);

    ifstream counts;
    counts.open(argv[1]);
    if (!counts) {
        cerr << "Unable to open the input file" << argv[1];
        return 1;
    }

    int total=0;
    map<string, int> pegcounts;
    int c; string peg;
    // read each line and split into two tokens
    while (counts >> c >> peg)
    {
        pegcounts[peg] = c;
        total += c;
    }
    counts.close();

    ostringstream p;
    p << argv[2] << '/' << ofname;
    ofstream outs(p.str());
    if (!outs) {
        cerr << "Unable to open the output file " << p.str() << endl;
        argv[1];
        return 2;
    }

    for( auto const& x : pegcounts )
    {
        outs << x.first << " " << (float) x.second/ (float) total << endl;
    }
    outs.close();

    return 0;
}

