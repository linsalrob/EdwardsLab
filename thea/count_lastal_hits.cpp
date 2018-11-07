//
// Created by redwards on 10/20/18.
// Process every file in the lastal directory and make a count of all the hits
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <map>
#include <unordered_set>
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

// FILTER THIS SO THAT WE DON'T COUNT THE SAME GENOME TWICE PER READ
// filter with 1e-5 evalue

    cout << "Filename\tpeg\tcounts\ttotal per mg\tnormalized" << endl;
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
        map<string, unordered_set<string>> seen;
        int total;
        while (getline(instream, line)) {
            if ('#' != line[0]) {
                istringstream ss(line);
                // Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, query length, subject length
                string peg; string seq_read; float per_ident; int aln_len; int m_m; int open_s; int q_start; int q_end; int s_start; int s_end;
                float e_val; int bit_score; int q_len; int s_len;
                ss >> seq_read >> peg >> per_ident >> aln_len >> m_m >> open_s >> q_start >> q_end >> s_start >> s_end >> e_val >> bit_score >> q_len >> s_len;

                // ignore everything with an e value above 1e-5
                if (1e-5 > e_val)
                    continue;

                // split the genome.orf to genome and orf
                vector<string> parts;
                boost::split(parts, peg, [](char c){return '.' == c;});

                // check to see if we have seen this read mapping to this genome before. If so, ignore it
                // If not, add the genome to the list per read.
                if (seen.count(seq_read) > 0 && seen[seq_read].count(parts[0]) > 0)
                    continue;
                if (0 == seen.count(seq_read))
                    seen[seq_read] = unordered_set<string>();
                seen[seq_read].insert(parts[0]);

                counts[peg]++;
                total++;
            }
        }


        boost::iostreams::close(inbuf);
        file.close();

        for (auto const& x : counts)
            cout << filename.path().string() << "\t" <<  x.first << "\t" << x.second << "\t" << total << "\t" << float(x.second) / total << endl;

    }
}






