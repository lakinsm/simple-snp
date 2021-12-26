#ifndef ASFFAST_PARSER_JOB_H
#define ASFFAST_PARSER_JOB_H

#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>
#include "concurrent_buffer_queue.h"
#include "args.h"


class ParserJob {
public:
    ParserJob(const std::string &parameter_string,
              const std::string &output_dir,
              ConcurrentBufferQueue* buffer_q,
              Args &args);
    ~ParserJob();

    void printInfo();
    void run();

    std::string sam_filepath;
    std::string samplename;
    std::string sam_readgroup;
    std::string sam_sampleid;
    std::string reference_name;
    std::string this_parent_ref;
    std::vector< std::string > this_children_ref;
    std::unordered_map< std::string, std::vector< std::vector< int > > > nucleotide_counts;
    std::unordered_map< std::string, std::vector< std::vector< long > > > qual_sums;
    std::unordered_map< std::string, std::vector< std::vector< long > > > mapq_sums;
    std::vector< long > ref_lens;

private:
    Args& _args;
    ConcurrentBufferQueue* _buffer_q;
    std::string _output_dir;

    void _addAlignedRead(const std::string &ref,
                         const std::string &cigar,
                         const std::string &seq,
                         const std::string &qual,
                         const long &pos,
                         const int &mapq);
    std::vector< std::string > _parseSamLine(const std::string &sam_line);
    void _writePositionalData();
    const std::unordered_map< char, int > _iupac_map = {
            {'A', 0},
            {'C', 1},
            {'G', 2},
            {'T', 3}
    };
};

#endif //ASFFAST_PARSER_JOB_H
