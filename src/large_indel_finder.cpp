#include "large_indel_finder.h"
#include "IntervalTree.h"
#include <limits>


typedef IntervalTree<long, int> ITree;


LargeIndelFinder::LargeIndelFinder(Args &args) : _args(args)
{

}


void LargeIndelFinder::findLargeIndels(const std::unordered_map< std::string,
                                       std::unordered_map< std::string,
                                       std::vector< std::vector< int > > > > &nucleotide_counts)
{
    // Find candidate ranges in each sample and order them by ascending size in a vector
    // Write this list out
    std::ofstream ofs1(_args.output_dir + "/large_indels.csv");
    ofs1 << "Sample,Reference,ReferenceAvgCoverage,Type,Start,Stop,RegionAvgCoverage,LeftBorderSharp,RightBorderSharp";
    ofs1 << std::endl;

    int range_idx = 0;
    std::vector< GenomicRange > all_ranges;
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
            std::vector< std::pair< long, long > > ref_ranges;
            std::vector< double > region_covs;
            std::vector< bool > range_high_confidence;
            _determineRanges(out_prefix, nucl, ofs1, ref_ranges, region_covs, range_high_confidence);
            for(int r = 0; r < ref_ranges.size(); ++r) {
                GenomicRange this_range(range_idx,
                                        range_idx,
                                        this_ref,
                                        ref_ranges[r].first,
                                        ref_ranges[r].second,
                                        this_range.stop - this_range.start + 1,
                                        range_high_confidence[r]);
                range_idx++;
                all_ranges.push_back(this_range);
            }
        }
    }
    std::sort(all_ranges.begin(), all_ranges.end());
    ofs1.close();

    // Tabling the below section for now, probably easier to do in R
//
//    // Pop the smallest range, insert other ranges into the interval tree, extract intersecting
//    // ranges, re-build the prio vector with updated ranges, and repeat until all ranges have been consumed.
//    std::unordered_map< int, bool > parent_range_confidence_map;
//    for(int i = 0; i < all_ranges.size(); ++i) {
//        parent_range_confidence_map[all_ranges[i].parent_id] = all_ranges[i].high_confidence;
//    }
//    while(!all_ranges.empty()) {
//        std::unordered_map< int, int > child_range_idx_map;
//        long query_start = all_ranges[0].start;
//        long query_stop = all_ranges[0].stop;
//        for(int i = 0; i < all_ranges.size(); ++i) {
//            child_range_idx_map[all_ranges[i].id] = i;
//        }
//        ITree::interval_vector interval_vec;
//        for(int i = 1; i < all_ranges.size(); ++i) {
//            interval_vec.push_back(ITree::interval(all_ranges[i].start, all_ranges[i].stop, all_ranges[i].id));
//        }
//        IntervalTree<long, int> interval_tree(std::move(interval_vec));
//
//        ITree::interval_vector results = interval_tree.findOverlapping(all_ranges[0].start, all_ranges[0].stop);
//        std::cout << "Query: " << all_ranges[0].start << ", " << all_ranges[0].stop << "\tSize: " << results.size() << std::endl;
//
//        // If the query is non-overlapping, extract and continue
//        if(results.size() == 0) {
//            // TODO
//            continue;
//        }
//
//        // Excise the query range, reformat the intersecting ranges, fill the range vector, and repeat
//        range_idx = 0;
//        std::vector< GenomicRange > new_all_ranges;
//
//        // First pass determines the nature of the overlaps detected by the interval tree
//        long from_right_min_idx = std::numeric_limits<long>max();
//        long from_left_max_idx = -1;
//        for(int i = 0; i < results.size(); ++i) {
//            if((results[i].start > query_start) && (results[i].start < from_right_min_idx)) {
//                from_right_min_idx = results[i].start;
//            }
//            if((results[i].stop < query_stop) && (results[i].stop > from_left_max_idx)) {
//                from_left_max_idx = results[i].stop;
//            }
//        }
//        if((from_left_max_idx < query_start) && (from_right_min_idx > query_stop)) {
//            // The query is entirely contained within the result intervals
//            std::string feature_name =
//        }
//        else {
//            if((from_left_max_idx > query_start) && (from_right_min_idx < query_stop)) {
//                // The query is partially overlapped from both sides: split into three queries and requeue
//                // TODO
//            }
//            else if((from_left_max_idx > query_start) && (from_right_min_idx > query_stop)) {
//                // The query is partially overlapped from the left only: split into two queries and requeue
//            }
//            else if((from_left_max_idx < query_start) && (from_right_min_idx < query_stop)) {
//                // The query is partially overlapped from the right only: split into two queries and requeue
//            }
//        }
//
//        for(int i = 0; i < results.size(); ++i) {
//            std::cout << '\t' << results[i].start << ", " << results[i].stop << ", " << results[i].value << std::endl;
//            bool left_overhang = results[i].start < all_ranges[0].start;
//            bool right_overhang = results[i].stop > all_ranges[0].stop;
//            if(left_overhang) {
//                GenomicRange this_range(all_ranges[child_range_idx_map[results[i].value]]);
//                this_range.id = range_idx++;
//                this_range.stop = results[i].start - 1;
//                this_range.size = this_range.stop - this_range.start + 1;
//                new_all_ranges.push_back(this_range);
//            }
//            if(right_overhang) {
//                GenomicRange this_range(all_ranges[child_range_idx_map[results[i].value]]);
//                this_range.id = range_idx++;
//                this_range.start = results[i].stop + 1;
//                this_range.size = this_range.stop - this_range.start + 1;
//                new_all_ranges.push_back(this_range);
//            }
//            if(left_overhang && right_overhang) {
//
//            }
//        }
//    }
}


void LargeIndelFinder::_determineRanges(const std::string &out_prefix,
                                        const std::vector< std::vector< int > > &nucl,
                                        std::ofstream &this_ofs,
                                        std::vector< std::pair< long, long > > &ranges,
                                        std::vector< double > &coverages,
                                        std::vector< bool > &high_confidence)
{
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
//                    std::cout << '\t' << j << '\t' << window_idx;
//                    std::cout << "\tloc: " << loc_bool_l << ',' << (!loc_bool_r);
//                    std::cout << " (" << this_depth << ", " << this_window_depth << ')';
//                    std::cout << "\twindow: " << window_bool_l << ',' << !(window_bool_r);
//                    std::cout << " (" << l_accel_avg << ", " << r_accel_avg << ')';
//                    std::cout << "\tborder: " << border_bool_l << ',' << border_bool_r << " (";
//                    std::cout << l_prev_ratio << ", " << r_prev_ratio << ')' << std::endl;
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
                    ranges.push_back(std::make_pair((long)j, (long)(j + window_idx)));
                    coverages.push_back(avg_region_depth);
                    if(border_bool_l && border_bool_r) {
                        high_confidence.push_back(true);
                    }
                    else {
                        high_confidence.push_back(false);
                    }
                }
            }
            j += window_idx;
        }
        else {
            prev_depth = this_depth;
        }
    }
}
#pragma clang diagnostic pop

#pragma clang diagnostic pop