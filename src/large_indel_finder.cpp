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
    // Write this list out
    std::ofstream ofs1(_args.output_dir + "/large_indels.csv");
    ofs1 << "Sample,Reference,ReferenceAvgCoverage,Type,Start,Stop,RegionAvgCoverage,LeftBorderSharp,RightBorderSharp";
    ofs1 << std::endl;
    for(auto &[sample, ref_map] : nucleotide_counts) {
        for(auto &[this_ref, nucl] : ref_map) {
            long total_ref_depth = 0;
            for(int i = 0; i < nucl.size(); ++i) {
                for(int j = 0; j < nucl[i].size(); ++j) {
                    total_ref_depth += nucl[i][j];
                }
            }
            double avg_ref_cov = (double)total_ref_depth / (double)nucl[0].size();
            std::string out_prefix = sample + ',' + this_ref + ',' + std::to_string(avg_ref_cov) + ',';
            std::cout << sample << '\t' << this_ref << std::endl;
            std::vector< std::pair< long, long > > this_ref_ranges;
            this_ref_ranges = _determineRanges(out_prefix, nucl, ofs1);
        }
    }
    ofs1.close();

    // Pop the smallest range, insert other ranges into the interval tree, extract intersecting
    // ranges, re-build the prio queue with updates ranges, and repeat until all ranges have been consumed.

}


std::vector< std::pair< long, long > > LargeIndelFinder::_determineRanges(std::string &out_prefix,
                                                                          const std::vector< std::vector< int > > &nucl,
                                                                          std::ofstream &this_ofs)
{
    std::vector< std::pair< long, long > > return_values;
    int ref_len = nucl[0].size();
    int prev_depth = 0;
    for(int j = 0; j < ref_len; ++j) {
        int this_depth = 0;
        for(int i = 0; i < nucl.size(); ++i) {
            this_depth += nucl[i][j];
        }
        int l_accel_depth = this_depth;
        double l_accel_avg;
        int l_accel_len = 1;
        for(int k = 1; k < _args.indel_accel_window_size; ++k) {
            if((j + k) < ref_len) {
                l_accel_len++;
                for(int i = 0; i < nucl.size(); ++i) {
                    l_accel_depth += nucl[i][j + k];
                }
            }
        }
        l_accel_avg = (double)l_accel_depth / (double)l_accel_len;
        double l_prev_ratio;
        if(this_depth != 0) {
            if(prev_depth != 0) {
                l_prev_ratio = (double)this_depth / (double)prev_depth;
            }
            else {
                l_prev_ratio = std::numeric_limits<double>::max();
            }
        }
        else {
            l_prev_ratio = 0;
        }
        bool loc_bool_l, loc_bool_r, window_bool_l, window_bool_r, border_bool_l, border_bool_r;
        border_bool_l = (prev_depth > 0) && (l_prev_ratio <= _args.large_indel_border_ratio);
        loc_bool_l = this_depth <= _args.large_indel_max_window_depth;
        window_bool_l = l_accel_avg <= _args.large_indel_max_window_depth;
        if((loc_bool_l && window_bool_l) || border_bool_l) {
            int r_accel_depth;
            double r_accel_avg;
            int this_window_depth = this_depth;
            int window_idx = 0;
            long total_depth = this_depth;
            double r_prev_ratio = 0;
            while(true) {
                if((j + window_idx + 1) == ref_len) {
                    break;
                }
                window_idx++;
                this_window_depth = 0;
                for(int i = 0; i < nucl.size(); ++i) {
                    this_window_depth += nucl[i][j + window_idx];
                }
                r_accel_depth = this_window_depth;
                total_depth += this_window_depth;
                int r_accel_len = 1;
                for(int k = 1; k < _args.indel_accel_window_size; ++k) {
                    if((j + window_idx - k) >= 0) {
                        r_accel_len++;
                        for(int i = 0; i < nucl.size(); ++i) {
                            r_accel_depth += nucl[i][j + window_idx - k];
                        }
                    }
                }
                r_accel_avg = (double)r_accel_depth / (double)r_accel_len;
                if(this_window_depth != 0) {
                    if(prev_depth != 0) {
                        r_prev_ratio = (double)prev_depth / (double)this_window_depth;
                    }
                    else {
                        r_prev_ratio = 0;
                    }
                }
                else {
                    r_prev_ratio = std::numeric_limits<double>::max();
                }
                loc_bool_r = this_window_depth <= _args.large_indel_max_window_depth;
                window_bool_r = r_accel_avg <= _args.large_indel_max_window_depth;
                border_bool_r = (this_window_depth > _args.large_indel_max_window_depth)
                        && (r_prev_ratio <= _args.large_indel_border_ratio);

                if(border_bool_r) {
                    total_depth -= this_window_depth;
                    window_idx--;
                    break;
                }
                if((!window_bool_r) && (!loc_bool_r)) {
                    total_depth -= this_window_depth;
                    window_idx--;
                    break;
                }
                prev_depth = this_window_depth;
            }
            prev_depth = this_window_depth;
            if(window_idx >= _args.min_large_indel_len) {
                double avg_region_depth = (double)total_depth / (double)window_idx;
                if(avg_region_depth < _args.large_indel_avg_max_depth) {
                    std::cout << '\t' << j << '\t' << window_idx;
                    std::cout << "\tloc: " << loc_bool_l << ',' << (!loc_bool_r);
                    std::cout << " (" << this_depth << ", " << this_window_depth << ')';
                    std::cout << "\twindow: " << window_bool_l << ',' << !(window_bool_r);
                    std::cout << " (" << l_accel_avg << ", " << r_accel_avg << ')';
                    std::cout << "\tborder: " << border_bool_l << ',' << border_bool_r << " (";
                    std::cout << l_prev_ratio << ", " << r_prev_ratio << ')' << std::endl;
                    this_ofs << out_prefix;
                    this_ofs << "deletion," << (j+1) << ',' << (j + window_idx + 1) << ',';
                    this_ofs << std::to_string(avg_region_depth) << ',';
                    if(border_bool_l) {
                        this_ofs << "TRUE,";
                    }
                    else {
                        this_ofs << "FALSE,";
                    }
                    if(border_bool_r) {
                        this_ofs << "TRUE";
                    }
                    else {
                        this_ofs << "FALSE";
                    }
                    this_ofs << std::endl;
                    return_values.push_back(std::make_pair((long)j, (long)(j + window_idx)));
                }
            }
            j += window_idx;
        }
        else {
            prev_depth = this_depth;
        }
    }
    return return_values;
}
