#ifndef SIMPLE_SNP_FASTA_PARSER_H
#define SIMPLE_SNP_FASTA_PARSER_H

#include <vector>
#include <string>
#include <unordered_map>


class FastaParser {
public:
    FastaParser(const std::string &fasta_filepath);

    void parseFasta(const std::vector< std::string > &selected_headers);

    std::string header;
    std::string seq;
    std::unordered_map< std::string, std::string > headers_seqs;
    std::unordered_map< std::string, long > headers_lens;
private:
    std::string _fasta_path;
};


#endif //SIMPLE_SNP_FASTA_PARSER_H
