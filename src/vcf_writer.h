#ifndef SIMPLE_SNP_VCF_WRITER_H
#define SIMPLE_SNP_VCF_WRITER_H
#include <chrono>
#include <ctime>
#include <string>
#include <vector>
#include <map>
#include <fstream>


typedef std::chrono::system_clock Clock;


struct vcfLineData {
    std::string chrom;
    std::string ref;
    std::vector< std::string > alt;
    long pos;
    double qual;
    int dp;
    int ns;
    int nsa;
    std::vector< int > alt_ns;
    std::vector< int > ac;
    std::vector< double > af;
    int ro;
    std::vector< int > ao;
    int ao_sum;
    std::vector< double > mqm;
    std::vector< double > gq;
    double mqmr;
    std::vector< std::string > type;
    std::vector< std::string > cigar;
};


class VcfWriter {
public:
    VcfWriter(std::string &vcf_path);

    void writeHeaders(const std::string &reference_path,
                      const std::string &commandline,
                      const std::vector< std::string > &contig_names,
                      const std::vector< long > &contig_lens);
    void writeSamples(const std::vector< std::string > &samplenames);
    void writeSampleData(const vcfLineData &vcf_line_data,
                         const std::map< std::string, std::string > &vcf_variants);
    void open();
    void close();

private:
    std::string _vcf_path;
    std::ofstream _ofs;
    std::vector< std::string > _sample_order;
};


#endif //SIMPLE_SNP_VCF_WRITER_H
