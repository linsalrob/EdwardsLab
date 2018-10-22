//
// Created by redwards on 10/20/18.
// Process every file in the lastal directory and make a count of all the hits
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <map>
#include <zlib.h>
using namespace std;

int main (int argc, char* argv[]) {

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <directory to process>\n";
        return 1;
    }

    if (!boost::filesystem::is_directory(argv[1])) {
        cerr << "Please provide a directory with all the files to summarize" << endl;
        return 1;
    }

    map<string, float> counts;
    map<string, int> total;
    for (auto &filename : boost::make_iterator_range(boost::filesystem::directory_iterator(argv[1]), {})) {

        ifstream file(filename.path().string(), ios_base::in | ios_base::binary);
        boost::iostreams::filtering_streambuf <boost::iostreams::input> inbuf;
        inbuf.push(boost::iostreams::gzip_decompressor());
        inbuf.push(file);

        //Convert streambuf to istream
        istream instream(&inbuf);


        string line;
        map<string, float> counts;
        int total;
        while (getline(instream, line)) {
            if ('#' != line[0]) {
                istringstream ss(line);
                string peg; string match;
                ss >> match >> peg;
                counts[peg]++;
                total++;
            }
        }


        boost::iostreams::close(inbuf);
        file.close();

        for (auto const& x : counts)
            cout << filename.path().string() << "\t" <<  x.first << "\t" << float(x.second) / total << endl;

    }
}






