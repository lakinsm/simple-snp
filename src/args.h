#ifndef ASFFAST_ARGS_H
#define ASFFAST_ARGS_H

#include <algorithm>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>
#include <limits.h>


class Args {
public:
    Args(int argc, const char *argv[]);

    std::string sam_file_dir;
    std::string output_dir;
    std::string reference_path;
    int min_intra_sample_alt;
    int min_inter_sample_alt;
    int min_intra_sample_depth;
    int min_inter_sample_depth;
    double min_major_freq;
    double min_minor_freq;
    int threads;

private:
    std::string _findFullDirPath(std::string path);
    void _usage();
};


#endif //ASFFAST_ARGS_H
