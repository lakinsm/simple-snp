#ifndef SIMPLE_SNP_LARGE_INDEL_FINDER_H
#define SIMPLE_SNP_LARGE_INDEL_FINDER_H

#include "args.h"
#include <vector>
#include <unordered_map>
#include <string>
#include <utility>
#include <algorithm>
#include <fstream>


class LargeIndelFinder {
public:
    LargeIndelFinder(Args &args);

    void findLargeIndels(const std::unordered_map< std::string,
                         std::unordered_map< std::string,
                         std::vector< std::vector< int > > > > &nucleotide_counts);

private:
    std::vector< std::pair< long, long > > _determineRanges(const std::vector< std::vector< int > > &nucl,
                                                            std::ofstream &this_ofs);

    Args& _args;
    std::unordered_map< std::string, std::pair< long, std::vector< std::string > > > _large_indels;
};


#endif //SIMPLE_SNP_LARGE_INDEL_FINDER_H
