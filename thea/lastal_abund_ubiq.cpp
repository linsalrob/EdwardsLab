/**
Created by redwards on 10/21/18.

This reads the output from count_lastal_hits.cpp which is a simple program
that generates a tuple of [mg id, peg, normalized abundance] and (a) averages the abundance
and (b) counts the metagenomes


 */


#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <map>
using namespace std;

int main (int argc, char* argv[]) {

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <file of count_lastal output>\n";
        return 1;
    }


    map<string, float> counts;
    map<string, int> total;

    ifstream reader;
    reader.open(argv[1]);
    if (!reader) {
        cerr << "Unable to open the input file" << argv[1] << endl;
        return 1;
    }

    int n; int tot; float c; string peg; string mg;
    // read each line and split into five tokens
    @TODO NEED TO CONVERT THIS TO IGNORE FIRST LINE!

    while (reader >> mg) {
        if ("Filename" == mg) {
            while ("normalized" != mg)
                reader >> mg;
            reader >> mg;
        }
        reader >> peg >> n >> tot >> c;
        counts[peg] += c;
        total[peg]++;
    }

    if (reader.fail())
        cerr << "Reading the file " << argv[1] << " failed" << endl;

    reader.close();

    for( auto const& x : counts )
    {
        cout << x.first << " " << (float) x.second/ (float) total[x.first] << " " << total[x.first] << endl;
    }

}