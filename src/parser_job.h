#ifndef ASFFAST_PARSER_JOB_H
#define ASFFAST_PARSER_JOB_H

#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>
#include "concurrent_buffer_queue.h"


class ParserJob {
public:
    ParserJob(const std::string &parameter_string, ConcurrentBufferQueue* buffer_q);
    ~ParserJob();

    void printInfo();
    void run();

    std::string sam_filepath;
    std::string samplename;
    std::string sam_readgroup;
    std::string sam_sampleid;
    std::string reference_name;
    std::vector< std::vector< int > > nucleotide_counts;
    std::vector< std::vector< long > > qual_sums;
    std::vector< std::vector< long > > mapq_sums;
    long ref_len;

private:
    ConcurrentBufferQueue* _buffer_q;

    void _addAlignedRead(const std::string &cigar,
                         const std::string &seq,
                         const std::string &qual,
                         const long &pos,
                         const int &mapq);
    std::vector< std::string > _parseSamLine(const std::string &sam_line);
    const std::unordered_map< char, int > _iupac_map = {
            {'A', 0},
            {'C', 1},
            {'G', 2},
            {'T', 3}
    };
};

#endif //ASFFAST_PARSER_JOB_H
