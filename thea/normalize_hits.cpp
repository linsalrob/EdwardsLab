//
// Created by redwards on 10/3/18.
//
// Normalize the rapsearch hits. This just sums the counts and then prints the numbers
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
using namespace std;

int main (int argc, char* argv[]) {

    if ( argc < 3) {
        cerr << "Usage: " << argv[0] << " <file to normalize> <directory to write to>\n";
        return 1;
    }


    ifstream counts;
    counts.open(argv[1]);
    if (!counts) {
        cerr << "Unable to open the input file" << argv[1];
        return 1;
    }

    string line;
    ostringstream p;
    p << argv[2] << '/' << argv[1] << endl;
    ofstream outs(p.str());
    if (!outs) {
        cerr << "Unable to open the output file " << p.str();
        argv[1];
        return 2;
    }


    int c=0;

    while (counts >> line)
    {
        // note this makes a stream from the string
        istringstream iss(line);
        cout << " read: " << line;
        int c;
        string peg;
        iss >> c;
        iss >> peg;
        cout << " peg : " << peg << " count: " << c;

    }

    return 0;
}

