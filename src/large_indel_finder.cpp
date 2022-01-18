#include "large_indel_finder.h"
#include "IntervalTree.h"

#include <limits>


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
    int ref_len = nucl[0].size();
    for(int j = 0; j < ref_len; ++j) {
        int this_depth = 0;
        int prev_depth = 0;
        for(int i = 0; i < nucl.size(); ++i) {
            this_depth += nucl[i][j];
        }
        double this_prev_ratio = 0;
        if(this_depth != 0) {
            this_prev_ratio = (double)this_depth / (double)prev_depth;
        }
//        std::cout << '\t' << j << '\t' << this_depth << std::endl;
        if((this_depth <= _args.large_indel_max_window_depth) ||
           ((this_prev_ratio > 0) && (this_prev_ratio <= _args.large_indel_border_ratio))) {
//            std::cout << '\t' << j << '\t' << "Window Trigger" << std::endl;
            long total_depth = (long)this_depth;
            int this_window_depth = this_depth;
            int window_idx = 0;
            while(true) {
                if((j + window_idx + 1) == ref_len) {
                    break;
                }
                window_idx++;
                this_window_depth = 0;
                for(int i = 0; i < nucl.size(); ++i) {
                    this_window_depth += nucl[i][j + window_idx];
                }
                total_depth += this_window_depth;
                if(this_window_depth != 0) {
                    if(prev_depth != 0) {
                        this_prev_ratio = (double)prev_depth / (double)this_window_depth;
                    }
                    else {
                        this_prev_ratio = std::numeric_limits<double>::max();
                    }
                }
                else {
                    this_prev_ratio = 0;
                }
                bool loc_bool = this_window_depth <= _args.large_indel_max_window_depth;
                bool window_bool = ((double)total_depth / (double)window_idx) <= _args.large_indel_max_window_depth;
                bool border_bool = (this_window_depth > 0) && (this_prev_ratio <= _args.large_indel_border_ratio);

                if(border_bool) {
                    std::cout << "\t\tWINDOW\t" << (j + window_idx) << "\tloc: " << loc_bool << " (" << this_window_depth << ')';
                    std::cout << "\twindow: " << window_bool << " (" << ((double)total_depth / (double)window_idx) << ')';
                    std::cout << "\tborder: " << border_bool << " (" << this_prev_ratio << ')' << std::endl;
                    window_idx--;
                    break;
                }
                if((!window_bool) && (!loc_bool)) {
                    std::cout << "\t\tWINDOW\t" << (j + window_idx) << "\tloc: " << loc_bool << " (" << this_window_depth << ')';
                    std::cout << "\twindow: " << window_bool << " (" << ((double)total_depth / (double)window_idx) << ')';
                    std::cout << "\tborder: " << border_bool << " (" << this_prev_ratio << ')' << std::endl;
                    window_idx--;
                    break;
                }
                prev_depth = this_window_depth;
            }
            prev_depth = this_window_depth;
//            std::cout << '\t' << j << '\t' << window_idx << '\t' << "CANDIDATE" << std::endl;
            if(window_idx >= _args.min_large_indel_len) {
                std::cout << '\t' << j << '\t' << window_idx << '\t' << "SELECTED" << std::endl;
                return_values.push_back(std::make_pair((long)j, (long)(j + window_idx)));
            }
            std::cout << "\t\t\tEnd idx: " << (j + window_idx) << std::endl;
            j += window_idx;  // TODO check
        }
        else {
            prev_depth = this_depth;
        }
    }
    return return_values;
}
