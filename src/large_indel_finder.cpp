#include "large_indel_finder.h"
#include "IntervalTree.h"


LargeIndelFinder::LargeIndelFinder(Args &args) : _args(args)
{

}


void LargeIndelFinder::findLargeIndels(const std::unordered_map< std::string,
                                       std::unordered_map< std::string,
                                       std::vector< std::vector< int > > > > &nucleotide_counts)
{
    // Find candidate ranges in each sample and order them by ascending size in a prio queue
    for(auto &[sample, ref_map] : nucleotide_counts) {
        for(auto &[this_ref, nucl] : ref_map) {
            std::cout << sample << '\t' << this_ref << std::endl;
            std::vector< std::pair< long, long > > this_ref_ranges;
            this_ref_ranges = _determineRanges(nucl);
        }
    }

    // Pop the smallest range, insert other ranges into the interval tree, extract intersecting
    // ranges, re-build the prio queue with updates ranges, and repeat until all ranges have been consumed.

}


std::vector< std::pair< long, long > > LargeIndelFinder::_determineRanges(const std::vector< std::vector< int > > &nucl)
{
    std::vector< std::pair< long, long > > return_values;
    for(int i = 0; i < nucl.size(); ++i) {
        int this_depth = 0;
        int prev_depth = 0;
        for(int j = 0; j < nucl[i].size(); ++j) {
            this_depth += nucl[i][j];
        }
        double this_prev_ratio = 0;
        if(this_depth != 0) {
            this_prev_ratio = (double)this_depth / (double)prev_depth;
        }
        if((this_depth <= _args.large_indel_max_window_depth) ||
           ((this_prev_ratio > 0) && (this_prev_ratio <= _args.large_indel_border_ratio))) {
            long total_depth = (long)this_depth;
            int this_window_depth = this_depth;
            int window_idx = 0;
            while(true) {
                if((i + window_idx + 1) == nucl.size()) {
                    break;
                }
                window_idx++;
                this_window_depth = 0;
                for(int j = 0; j < nucl[i].size(); ++j) {
                    this_window_depth += nucl[i + window_idx][j];
                }
                total_depth += this_window_depth;
                if(this_window_depth != 0) {
                    this_prev_ratio = (double)prev_depth / (double)this_window_depth;
                }
                else {
                    this_prev_ratio = 0;
                }
                bool loc_bool = this_window_depth <= _args.large_indel_border_ratio;
                bool window_bool = ((double)total_depth / (double)window_idx) <= _args.large_indel_border_ratio;
                bool border_bool = (this_prev_ratio > 0) && (this_prev_ratio <= _args.large_indel_border_ratio);
                if(border_bool) {
                    window_idx--;
                    break;
                }
                if((!window_bool) || (!loc_bool)) {
                    window_idx--;
                    break;
                }
                prev_depth = this_window_depth;
            }
            prev_depth = this_window_depth;
            std::cout << '\t' << i << '\t' << window_idx << '\t' << "CANDIDATE" << std::endl;
            if(window_idx >= _args.min_large_indel_len) {
                std::cout << '\t' << i << '\t' << window_idx << '\t' << "SELECTED" << std::endl;
                return_values.push_back(std::make_pair((long)i, (long)(i + window_idx));
            }
            std::cout << "\t\t\tEnd idx: " << (i + window_idx) << std::endl;
            i += window_idx;  // TODO check
        }
        else {
            prev_depth = this_depth;
        }
    }
    return return_values;
}
