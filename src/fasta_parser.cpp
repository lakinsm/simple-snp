#include "fasta_parser.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <set>


FastaParser::FastaParser(const std::string &fasta_filepath) : _fasta_path(fasta_filepath)
{

}


void FastaParser::parseFasta(const std::vector< std::string > &selected_headers)
{
    std::ifstream ifs;
    ifs.open(_fasta_path, std::ifstream::in);

    if(!ifs.is_open()) {
        std::cerr << "ERROR: Could not open FASTA file " << _fasta_path << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::set< std::string > selects(selected_headers.begin(), selected_headers.end());

    std::string temp_line;
    std::getline(ifs, temp_line);
    while(ifs.good()) {
        if(temp_line[0] != '>') {
            std::cerr << "ERROR: FASTA file headers must begin with >, provided: " << temp_line << std::endl;
            std::exit(EXIT_FAILURE);
        }
        header = "";
        seq = "";

        header = temp_line.substr(1, temp_line.find(' ') - 1);
        std::getline(ifs, temp_line);
        while((!ifs.eof()) && (temp_line[0] != '>')) {
            seq += temp_line;
            std::getline(ifs, temp_line);
        }
        if(headers_seqs.count(header)) {
            std::cerr << "ERROR: FASTA headers must be unqiue, duplicated provided: " << header << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if((!select.empty()) and (!selects.count(header))) {
            continue;
        }
        headers_seqs[header] = seq;
        headers_lens[header] = (long)seq.length();
    }

    ifs.close();
}
