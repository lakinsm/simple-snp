#ifndef SIMPLE_SNP_LARGE_INDEL_FINDER_H
#define SIMPLE_SNP_LARGE_INDEL_FINDER_H

#include "args.h"
#include <vector>
#include <unordered_map>
#include <string>
#include <utility>
#include <algorithm>
#include <fstream>


struct GenomicRange {
    GenomicRange(int this_id,
                 std::string this_ref,
                 long this_start,
                 long this_stop,
                 long this_size,
                 bool this_confidence)
                 : id(this_id),
                 ref(this_ref),
                 start(this_start),
                 stop(this_stop),
                 size(this_size),
                 high_confidence(this_confidence) {};

    bool operator < (const GenomicRange& other) const
    {
        return (size < other.size);
    }

    int id;
    std::string ref;
    long start;
    long stop;
    long size;
    bool high_confidence;
};


class LargeIndelFinder {
public:
    LargeIndelFinder(Args &args);

    void findLargeIndels(const std::unordered_map< std::string,
                         std::unordered_map< std::string,
                         std::vector< std::vector< int > > > > &nucleotide_counts);

private:
    void _determineRanges(const std::string &out_prefix,
                          const std::vector< std::vector< int > > &nucl,
                          std::ofstream &this_ofs,
                          std::vector< std::pair< long, long > > &ranges,
                          std::vector< double > &coverages,
                          std::vector< bool > &high_confidence);

    Args& _args;
    std::unordered_map< std::string, std::pair< long, std::vector< std::string > > > _large_indels;
};


#endif //SIMPLE_SNP_LARGE_INDEL_FINDER_H
