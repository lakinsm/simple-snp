#ifndef ASFFAST_PARSER_JOB_H
#define ASFFAST_PARSER_JOB_H

#include <string>
#include <iostream>
#include <vector>
#include <set>
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
    std::vector< std::string > contents;
    std::set< std::string > seen_headers;
    std::set< std::string > aligned_headers;
    long reads_processed;
    long reads_aligned;

private:
    ConcurrentBufferQueue* _buffer_q;

    std::vector< std::string > _parseSamLine(const std::string &sam_line);
};

#endif //ASFFAST_PARSER_JOB_H
