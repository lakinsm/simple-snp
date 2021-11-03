#ifndef SIMPLE_SNP_VCF_WRITER_H
#define SIMPLE_SNP_VCF_WRITER_H
#include <chrono>
#include <ctime>
#include <string>
#include <vector>


typedef std::chrono::system_clock Clock;


struct vcfLineData {
    std::string chrom;
    std::string ref;
    std::vector< std::string > alt;
    long pos;
    double qual;
    int dp;
    int ns;
    std::vector< int > alt_ns;
    int ac;
    std::vector< double > af;
    int ro;
    std::vector< int > ao;
    std::vector< double > mqm;
    double mqmr;
    std::vector< std::string > type;
    std::vector< std::string > cigar;
};


class VcfWriter {
public:
    VcfWriter(std::ofstream &ofs);

    void writeHeaders(const std::string &reference_path);
    void writeSamples(const std::vector< std::string > &samplenames);
    void writeSampleData(const vcfLineData &vcf_line_data,
                         const std::map< std::string, std::string > &vcf_variants);

private:
    std::ofstream _ofs;
    std::vector< std::string > _sample_order;
};


#endif //SIMPLE_SNP_VCF_WRITER_H
