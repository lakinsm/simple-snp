#include "vcf_writer.h"
#include <iostream>
#include <cassert>


VcfWriter::VcfWriter(std::string &vcf_path) : _vcf_path(vcf_path)
{

}


void VcfWriter::writeHeaders(const std::string &reference_path,
                             const std::string &commandline,
                             const std::vector< std::string > &contig_names,
                             const std::vector< long > &contig_lens)
{
    auto now = Clock::now();
    std::time_t now_c = Clock::to_time_t(now);
    struct tm *parts = std::localtime(&now_c);
    std::string time_string = "";
    time_string += std::to_string(1900 + parts->tm_year) + std::to_string(1 + parts->tm_mon);
    if(parts->tm_mday < 10) {
        time_string += "0" + std::to_string(parts->tm_mday);
    }
    else {
        time_string += std::to_string(parts->tm_mday);
    }

    _ofs << "##filterformat=VCFv4.2" << std::endl;
    _ofs << "##fileDate=" << time_string << std::endl;
    _ofs << "##source=SimpleSNP v0.1" << std::endl;
    _ofs << "##reference=" << reference_path << std::endl;
    for(int i = 0; i < contig_names.size(); ++i) {
        _ofs << "##contig=<ID=" << contig_names[i] << ",length=" << contig_lens[i] << ">" << std::endl;
    }
    _ofs << "##phasing=none" << std::endl;
    _ofs << "##commandline=" << commandline << std::endl;
    _ofs << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">" << std::endl;
    _ofs << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">" << std::endl;
    _ofs << "##INFO=<ID=NSA,Number=1,Type=Integer,Description=\"Number of samples with alternate alleles\">" << std::endl;
    _ofs << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">" << std::endl;
    _ofs << "##INFO=<ID=RO,Number=1,Type=Integer,Description=\"Total reference alleles observed\">" << std::endl;
    _ofs << "##INFO=<ID=AO,Number=A,Type=Integer,Description=\"Total alternate alleles observed\">" << std::endl;
    _ofs << "##INFO=<ID=MQM,Number=A,Type=Float,Description=\"Mean mapping quality of alternate alleles observed\">";
    _ofs << std::endl;
    _ofs << "##INFO=<ID=MQMR,Number=1,Type=Float,Description=\"Mean mapping quality of reference allele\">" << std::endl;
    _ofs << "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"Type of allele, either snp, mnp, ins, del, or copmlex\">";
    _ofs << std::endl;
    _ofs << "##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"The extended CIGAR representation of each ";
    _ofs << "alternate allele\">" << std::endl;
    _ofs << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\"" << std::endl;
    _ofs << "##FORMAT=<ID=GQ,Number=A,Type=Float,Description=\"Genotype quality (avgerage PHRED score of alleles observed)\"";
    _ofs << std::endl;
    _ofs << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\"" << std::endl;
    _ofs << "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Number of observation for each allele\"" << std::endl;
    _ofs << "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype likelihood\"" << std::endl;
    _ofs << "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observation count\"" << std::endl;
    _ofs << "##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observation count\"" << std::endl;
    _ofs << "##FORMAT=<ID=QR,Number=1,Type=Float,Description=\"Sum of reference allele PHRED scores\"" << std::endl;
    _ofs << "##FORMAT=<ID=QA,Number=A,Type=Float,Description=\"Sum of alternate allele PHRED scores\"" << std::endl;
}


void VcfWriter::writeSamples(const std::vector< std::string > &samplenames)
{
    _sample_order = samplenames;
    _ofs << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for(int i = 0; i < samplenames.size(); ++i) {
        _ofs << '\t' << samplenames[i];
    }
    _ofs << std::endl;
}


void VcfWriter::writeSampleData(const vcfLineData &vcf_line_data,
                                const std::map< std::string, std::string > &vcf_variants)
{
    _ofs << vcf_line_data.chrom;
    _ofs << '\t' << std::to_string(vcf_line_data.pos);
    _ofs << "\t.\t" << vcf_line_data.ref;
    _ofs << '\t' << vcf_line_data.alt[0];
    for(int i = 1; i < vcf_line_data.alt.size(); ++i) {
        _ofs << ',' << vcf_line_data.alt[i];
    }
    _ofs << '\t' << std::to_string(vcf_line_data.qual) << "\t.\t";
    _ofs << "NSA=" << std::to_string(vcf_line_data.nsa) << ';';
    _ofs << "AC=" << std::to_string(vcf_line_data.ac) << ';';
    _ofs << "AF=" << std::to_string(vcf_line_data.af[0]);
    for(int i = 1; i < vcf_line_data.af.size(); ++i) {
        _ofs << ',' << std::to_string(vcf_line_data.af[i]);
    }
    _ofs << ';';
    _ofs << "AO=" << std::to_string(vcf_line_data.ao[0]);
    _ofs << "RO=" << std::to_string(vcf_line_data.ro) << ';';
    for(int i = 1; i < vcf_line_data.ao.size(); ++i) {
        _ofs << ',' << std::to_string(vcf_line_data.ao[i]);
    }
    _ofs << ';';
    _ofs << "CIGAR=1X;";
    _ofs << "DP=" << std::to_string(vcf_line_data.dp) << ';';
    _ofs << "MQM=" << std::to_string(vcf_line_data.mqm[0]);
    for(int i = 1; i < vcf_line_data.mqm.size(); ++i) {
        _ofs << ',' << std::to_string(vcf_line_data.mqm[i]);
    }
    _ofs << ';';
    _ofs << "MQMR=" << std::to_string(vcf_line_data.mqmr) << ';';
    _ofs << "NS=" << std::to_string(_sample_order.size()) << ';';
    _ofs << "TYPE=snp\t";
    _ofs << "GT:DP:AD:RO:QR:AO:QA:GL";
    for(int i = 0; i < _sample_order.size(); ++i) {
        _ofs << '\t' << vcf_variants.at(_sample_order[i]);
    }
    _ofs << std::endl;
}

void VcfWriter::open()
{
    _ofs.open(_vcf_path);
    if(!_ofs.is_open()) {
        std::cerr << "ERROR: VCF ofstream handle broken." << std::endl;
        std::exit(EXIT_FAILURE);
    }
}


void VcfWriter::close()
{
    _ofs.close();
}
