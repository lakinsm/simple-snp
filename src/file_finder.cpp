#include "file_finder.h"
#include <glob.h>


std::vector< std::string > FileFinder::findSamFiles(const std::string &input_path)
{
    std::vector< std::string > return_files;
    glob_t glob_result;
    memset(&glob_result, 0, sizeof(glob_result));

    std::string glob_pattern = input_path + "/*.sam";
    int return_value = glob(glob_pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    if(return_value != 0 && return_value != 3) {
        globfree(&glob_result);
        std::cerr << "findSamFiles() glob() failed with return value: " << return_value << std::endl;
        std::exit(EXIT_FAILURE);
    }
    else if(return_value == 3) {
        std::cerr << std::endl << "The specified input directory (" << _args.test_input_directory;
        std::cerr << ") contains no detectable SAM files with extension .sam" << std::endl << std::endl;
        std::exit(EXIT_FAILURE);
    }
    else if(return_value == 0) {
        for(size_t i = 0; i < glob_result.gl_pathc; ++i) {
            std::string file(glob_result.gl_pathv[i]);
            return_files.push_back(file);
        }
    }
    return return_files;
}
