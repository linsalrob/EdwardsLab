//
// Created by redwards on 10/3/18.
// Process every file in the directory and make a count of all the hits
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <map>
using namespace std;

int main (int argc, char* argv[]) {

    if ( argc < 2) {
        cerr << "Usage: " << argv[0] << " <directory to process>\n";
        return 1;
    }

    if (! boost::filesystem::is_directory(argv[1])) {
        cerr << "Please provide a directory with all the files to summarize" << endl;
        return 1;
    }

    map<string, float> counts;
    map<string, int> total;
    for(auto& filename : boost::make_iterator_range(boost::filesystem::directory_iterator(argv[1]), {})) {
        ifstream reader;
        reader.open(filename.path().string());
        if (!reader) {
            cerr << "Unable to open the input file" << filename << endl;
            return 1;
        }
        float c; string peg;
        // read each line and split into two tokens
        while (reader >> peg >> c) {
            counts[peg] += c;
            total[peg]++;
        }
        reader.close();
    }

    for( auto const& x : counts )
    {
        cout << x.first << " " << (float) x.second/ (float) total[x.first] << " " << total[x.first] << endl;
    }

}






