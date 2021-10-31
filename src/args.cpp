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

    sam_file_dir = _findFullDirPath(arg_list[1]);
    output_dir = _findFullDirPath(arg_list[2]);
    threads = std::stoi(arg_list[3].c_str());

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
    std::cout << "\nUsage:" << std::endl;
    std::cout << "\tsam_parse_merge sam_file_dir/ output_dir/ threads";
    std::cout << std::endl << std::endl;
    std::exit(EXIT_FAILURE);
}
