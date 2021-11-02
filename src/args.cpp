#include <libgen.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <cassert>
#include "args.h"


// Public member functions
Args::Args(int argc, const char *argv[])
{
    std::vector< std::string > arg_list(argv, argv+argc);
    if(argc < 4) {
        std::cerr << std::endl << "Too few arguments." << std::endl;
        _usage();
    }

    // Default values
    min_intra_sample_alt = 3;
    min_intra_sample_depth = 5;
    min_inter_sample_alt = 7;
    min_inter_sample_depth = 10;
    min_major_freq = 0.5;
    min_minor_freq = 0.3;
    threads = 2;

    // User-specified required values
    sam_file_dir = _findFullDirPath(arg_list[1]);
    output_dir = _findFullDirPath(arg_list[2]);
    reference_path = arg_list[3];

    // User-specified optional values
    for(int i = 4; i < argc; ++i) {
        if(arg_list[i] == "-a")
            min_intra_sample_alt = std::stoi(arg_list[++i].c_str());
        else if(arg_list[i] == "-A")
            min_inter_sample_alt = std::stoi(arg_list[++i].c_str());
        else if(arg_list[i] == "-d")
            min_intra_sample_depth = std::stoi(arg_list[++i].c_str());
        else if(arg_list[i] == "-D")
            min_inter_sample_depth = std::stoi(arg_list[++i].c_str());
        else if(arg_list[i] == "-F")
            min_major_freq = std::stod(arg_list[++i].c_str());
        else if(arg_list[i] == "-f")
            min_minor_freq = std::stod(arg_list[++i].c_str());
        else if(arg_list[i] == "-t")
            threads = std::stoi(arg_list[++i].c_str());
        else {
            std::cerr << "ERROR: Option not recognized, " << arg_list[i] << std::endl;
            _usage();
        }
    }

    if(threads < 2) {
        std::cerr << "ERROR: Threads must be at least 2, provided: " << threads << std::endl;
        std::exit(EXIT_FAILURE);
    }
}


// Private member functions
std::string Args::_findFullDirPath(std::string path)
{
    char* symlinkpath = &path[0];
    char fullpath[PATH_MAX];
    char* ptr = realpath(symlinkpath, fullpath);
    std::string ret(ptr);
    return ret;
}


void Args::_usage()
{
    std::cout << std::endl << "Usage:" << std::endl;
    std::cout << "\tsimple_snp sam_file_dir/ output_dir/ reference.fasta [options]" << std::endl << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "\t-a\tWithin-sample minimum alternate allele count to call a variant [3]" << std::endl;
    std::cout << "\t-A\tBetween-sample minimum alternate allele count to call a variant [7]" << std::endl;
    std::cout << "\t-d\tWithin-sample minimum read depth to call a variant [5]" << std::endl;
    std::cout << "\t-D\tBetween-sample minimum read depth to call a variant [10]" << std::endl;
    std::cout << "\t-f\tMinimum within-sample alternate allele frequency to call a minor variant [0.3]" << std::endl;
    std::cout << "\t-f\tMinimum within-sample alternate allele frequency to call a major variant [0.5]" << std::endl;
    std::cout << "\t-t\tThreads to use, minimum 2 [2]" << std::endl;
    std::cout << std::endl << std::endl;
    std::exit(EXIT_FAILURE);
}
