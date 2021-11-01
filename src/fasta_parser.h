#ifndef SIMPLE_SNP_FASTA_PARSER_H
#define SIMPLE_SNP_FASTA_PARSER_H

#include <string>


class FastaParser {
public:
    FastaParser(const std::string &fasta_filepath);

    void parseFasta();

    std::string header;
    std::string seq;
private:
    std::string _fasta_path;
};


#endif //SIMPLE_SNP_FASTA_PARSER_H
