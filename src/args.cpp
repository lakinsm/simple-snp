#include <libgen.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <cassert>
#include <filesystem>
#include "args.h"


// Public member functions
Args::Args(int argc, const char *argv[])
{
    std::vector< std::string > arg_list(argv, argv+argc);
    if(argc < 4) {
        std::cerr << std::endl << "Too few arguments." << std::endl;
        _usage();
    }

    // User-specified required values
    sam_file_dir = _findFullDirPath(arg_list[1]);
    output_dir = arg_list[2];
    reference_path = _findFullDirPath(arg_list[3]);

    if(!std::filesystem::is_directory(output_dir)) {
        std::filesystem::create_directory(output_dir);
    }
    output_dir = _findFullDirPath(arg_list[2]);

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
        else if(arg_list[i] == "-n") {
            db_ann_file = _findFullDirPath(arg_list[++i]);
            std::size_t start_pos = db_ann_file.find(".ann");
            if(start_pos == std::string::npos) {
                std::cerr << "ERROR: Database annotation file (-n) must have .ann extension, provided: ";
                std::cerr << db_ann_file << std::endl;
                exit(EXIT_FAILURE);
            }
            db_names_file = db_ann_file;
            db_names_file.replace(start_pos, 4, ".names");
            if(!std::filesystem::exists(db_names_file)) {
                db_names_file = "";
            }
        }
        else {
            std::cerr << "ERROR: Option not recognized, " << arg_list[i] << std::endl;
            _usage();
        }
    }

    if(threads < 3) {
        std::cerr << "ERROR: Threads must be at least 3, provided: " << threads << std::endl;
        std::exit(EXIT_FAILURE);
    }
}


// Private member functions
std::string Args::_findFullDirPath(std::string path)
{
    if(!std::filesystem::exists(path)) {
        std::cerr << "ERROR: File does not exist: " << path << std::endl;
        exit(EXIT_FAILURE);
    }

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
    std::cout << "\t-n\tFILE\tComma-separated file linking reference ID to subregions of interest (.ann extension)";
    std::cout << std::endl;
    std::cout << "\t-t\tThreads to use, minimum 3 [3]" << std::endl;
    std::cout << std::endl << std::endl;
    std::exit(EXIT_FAILURE);
}
