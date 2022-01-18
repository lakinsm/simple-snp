#ifndef ASFFAST_ARGS_H
#define ASFFAST_ARGS_H

#include <algorithm>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>
#include <limits.h>
#include <unordered_map>


class Args {
public:
    Args(int argc, const char *argv[]);

    std::string sam_file_dir;
    std::string output_dir;
    std::string reference_path;
    std::string db_ann_file = "";
    std::string db_names_file = "";
    int min_intra_sample_alt = 3;
    int min_inter_sample_alt = 7;
    int min_intra_sample_depth = 5;
    int min_inter_sample_depth = 10;
    int min_large_indel_len = 48;
    int large_indel_max_window_depth = 2;
    int indel_accel_window_size = 3;
    double large_indel_border_ratio = 0.1;
    double min_major_freq = 0.7;
    double min_minor_freq = 0.4;
    int threads = 3;

    // { acc: < < start, stop, strand, gene, product > > }
    std::unordered_map< std::string, std::vector< std::vector< std::string > > > db_ann_map;

    // { acc_parent: < acc_child1, acc_child2, ... > }
    std::unordered_map< std::string, std::vector< std::string > > db_parent_map;

    // { acc_child: acc_parent }
    std::unordered_map< std::string, std::string > rev_db_parent_map;

    // { acc_parent: name }
    std::unordered_map< std::string, std::string > db_parent_name_map;

    // { acc_child: name }
    std::unordered_map< std::string, std::string > db_child_name_map;

private:
    std::string _findFullDirPath(std::string path);
    void _usage();
};


#endif //ASFFAST_ARGS_H
