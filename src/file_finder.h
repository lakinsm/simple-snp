#ifndef SIMPLE_SNP_FILE_FINDER_H
#define SIMPLE_SNP_FILE_FINDER_H

#include <string>
#include <vector>

class FileFinder {
public:
    FileFinder();

    std::vector< std::string > findSamFiles(const std::string &input_path);
};


#endif //SIMPLE_SNP_FILE_FINDER_H
