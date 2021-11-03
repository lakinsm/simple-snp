#include "fasta_parser.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <algorithm>


FastaParser::FastaParser(const std::string &fasta_filepath) : _fasta_path(fasta_filepath)
{

}


void FastaParser::parseFasta()
{
    std::ifstream ifs;
    ifs.open(_fasta_path, std::ifstream::in);

    if(!ifs.is_open()) {
        std::cerr << "Could not open FASTA file " << _fasta_path << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::string temp_line;
    std::getline(ifs, temp_line);
    if(temp_line[0] != '>') {
        std::cerr << "FASTA file headers must begin with >, provided: " << temp_line << std::endl;
        std::exit(EXIT_FAILURE);
    }

    header = temp_line.substr(1, temp_line.find(' '));
    std::getline(ifs, temp_line);
    while((!ifs.eof()) && (temp_line[0] != '>')) {
        seq += temp_line;
        std::getline(ifs, temp_line);
    }

    ifs.close();
}
