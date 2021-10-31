#include "parser_job.h"
#include <fstream>
#include <sstream>
#include <cassert>
#include <ctype.h>


ParserJob::ParserJob(const std::string &parameter_string,
                     const std::string &output_dir,
                     ConcurrentBufferQueue* buffer_q)
                     : _buffer_q(buffer_q), _output_dir(output_dir)
{
    std::stringstream ss;
    ss.str(parameter_string);
    std::getline(ss, sam_filepath, '|');
    std::getline(ss, samplename);
    ref_len = 0;
}


ParserJob::~ParserJob()
{
    _buffer_q->num_active_jobs -= 1;
    _buffer_q->num_completed_jobs += 1;
}


void ParserJob::printInfo()
{
    std::cout << std::endl;
    std::cout << samplename << '\t' << sam_sampleid << '\t' << sam_readgroup << '\t' << sam_filepath << std::endl;
    std::cout << reference_name << '\t' << ref_len << std::endl;
    std::cout << std::endl;
}


void ParserJob::run()
{
    std::string this_header, line;
    std::ifstream ifs(sam_filepath, std::ios::in);

    if(!ifs.good()) {
        return;
    }

    this_header = "";
    bool headers = false;
    bool readgroup_present = false;
    bool ref_info_present = false;
    while(!headers) {
        std::getline(ifs, line);

        if(line.empty()) {
            return;
        }

        if(line[0] == '@') {
            if(line.substr(0, 3) == "@SQ") {
                if(ref_info_present) {
                    std::cerr << "ERROR: Multiple reference contigs are not currently supported. More than one";
                    std::cout << " reference contig was found in SAM file: " << sam_filepath << std::endl;
                    std::exit(EXIT_FAILURE);
                }

                std::stringstream ss_sq;
                std::string sq_part;
                ss_sq.str(line);
                std::getline(ss_sq, sq_part, '\t');
                std::getline(ss_sq, sq_part, '\t');
                reference_name = sq_part.substr(3);
                std::getline(ss_sq, sq_part);
                std::string this_len_part = sq_part.substr(3);
                ref_len = std::stol(this_len_part.c_str());
                ref_info_present = true;
            }

            if(line.substr(0, 3) == "@RG") {
                std::stringstream ss_rg;
                std::string rg_part;
                ss_rg.str(line);
                std::getline(ss_rg, rg_part, '\t');
                std::getline(ss_rg, rg_part, '\t');
                sam_readgroup = rg_part.substr(3);
                std::getline(ss_rg, rg_part);
                sam_sampleid = rg_part.substr(3);
                readgroup_present = true;
            }
        }
        else {
            headers = true;
        }
    }

    if(!readgroup_present) {
        std::cerr << "ERROR: Readgroup information (@RG) is not present in SAM file (" << sam_filepath << ").";
        std::cerr << " @RG ID and SM must be set." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if(!ref_info_present) {
        std::cerr << "ERROR: Reference sequence information (@SQ) is not present in SAM file (" << sam_filepath << ").";
        std::cerr << " @SQ SN and LN must be present only once (one contig in reference file)." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if(ref_len <= 0) {
        std::cerr << "ERROR: Reference sequence length is not positive, provided: " << ref_len << std::endl;
        std::exit(EXIT_FAILURE);
    }

    for(int i = 0; i < _iupac_map.size(); ++i) {
        nucleotide_counts.push_back(std::vector< int >(ref_len, 0));
        qual_sums.push_back(std::vector< long >(ref_len, 0));
        mapq_sums.push_back(std::vector< long >(ref_len, 0));
    }

    std::vector< std::string > res;
    int sam_flag;
    res = _parseSamLine(line);
    if((res.size() == 0) || (res[0].empty())) {
        return;
    }
    sam_flag = std::stoi(res[0].c_str());
    if(((sam_flag & 4) == 0) and ((sam_flag & 256) == 0) and ((sam_flag & 2048) == 0)) {
        // Primary alignment
        _addAlignedRead(res[3], res[4], res[5], std::stol(res[1].c_str()), std::stoi(res[2].c_str()));
    }

    while(std::getline(ifs, line)) {
        res = _parseSamLine(line);
        sam_flag = std::stoi(res[0].c_str());
        if(((sam_flag & 4) == 0) and ((sam_flag & 256) == 0) and ((sam_flag & 2048) == 0)) {
            // Primary alignment
            _addAlignedRead(res[3], res[4], res[5], std::stol(res[1].c_str()), std::stoi(res[2].c_str()));
        }
    }

//    printInfo();

//    while(!_buffer_q->tryPush(contents, barcode, reads_processed, reads_aligned)) {}
}


void ParserJob::_addAlignedRead(const std::string &cigar,
                                const std::string &seq,
                                const std::string &qual,
                                const long &pos,
                                const int &mapq)
{
    long read_idx = 0;
    long target_idx = pos - 1;
    long cigar_idx = 0;
    std::string num = "";
    std::string op = "";
    while(cigar_idx != cigar.size()) {
        if(std::isdigit(cigar[cigar_idx])) {
            num += cigar[cigar_idx];
        }
        else {
            op = cigar[cigar_idx];
            int numeric_num = std::stoi(num.c_str());
            if((op == "M") or (op == "=") or (op == "X")) {
                for(int i = 0; i < numeric_num; ++i) {
                    nucleotide_counts[_iupac_map.at(seq[read_idx])][target_idx]++;
                    qual_sums[_iupac_map.at(seq[read_idx])][target_idx] += int(qual[read_idx]) - 33;  // Phred 33
                    mapq_sums[_iupac_map.at(seq[read_idx])][target_idx] += mapq;
                    read_idx++;
                    target_idx++;
                }
            }
            else if((op == "D") or (op == "N")) {
                target_idx += numeric_num;
            }
            else if((op == "I") or (op == "S")) {
                read_idx += numeric_num;
            }
            num = "";
        }
        cigar_idx++;
    }
}


std::vector< std::string > ParserJob::_parseSamLine(const std::string &sam_line)
{
    //      0             1           2     3     4    5
    // < sam flag, start pos 1-idx, mapq, cigar, seq, qual >
    std::vector< std::string > ret;
    std::stringstream this_ss;
    this_ss.str(sam_line);
    std::string this_entry;
    std::getline(this_ss, this_entry, '\t');  // 0. read name
    std::getline(this_ss, this_entry, '\t');  // 1. sam flag
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');  // 2. ref name
    std::getline(this_ss, this_entry, '\t');  // 3. start pos 1-idx
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');  // 4. mapq
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');  // 5. cigar
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');  // 6. rnext
    std::getline(this_ss, this_entry, '\t');  // 7. pnext
    std::getline(this_ss, this_entry, '\t');  // 8. tlen
    std::getline(this_ss, this_entry, '\t');  // 9. seq
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');  // 10. qual
    ret.push_back(this_entry);
    return ret;
}


void ParserJob::_writePositionalData()
{
    std::string outfile_path = _output_dir + "/" + samplename + "_positional_data.tsv";
    std::ofstream ofs(outfile_path);

    ofs << "ReferenceIndex\tA_count,A_avg_qual,A_avg_mapq\tC_count,C_avg_qual,C_avg_mapq\t";
    ofs << "G_count,G_avg_qual_G_avg_mapq\tT_count,T_avg_qual,T_avg_mapq" << std::endl;

    for(int j = 0; j < ref_len; ++j) {
        ofs << (i + 1);
        for(int i = 1; i < _iupac_map.size(); ++i) {
            ofs << "\t" << nucleotide_counts[i][j] << ",";
            ofs << ((double)qual_sums[i][j] / (double)nucleotide_counts[i][j]) << ",";
            ofs << ((double)mapq_sums[i][j] / (double)nucleotide_counts[i][j]);
        }
        ofs << std::endl;
    }
    ofs.close();
}
