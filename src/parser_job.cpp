#include "parser_job.h"
#include <fstream>
#include <sstream>
#include <cassert>
#include <ctype.h>


ParserJob::ParserJob(const std::string &parameter_string,
                     const std::string &output_dir,
                     ConcurrentBufferQueue* buffer_q,
                     Args &args)
                     : _buffer_q(buffer_q), _output_dir(output_dir), _args(args)
{
    std::stringstream ss;
    ss.str(parameter_string);
    std::getline(ss, sam_filepath, '|');
    std::getline(ss, samplename);
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
    for(int i = 0; i < this_children_ref.size(); ++i) {
        std::cout << '\t' << this_children_ref[i] << '\t' << ref_lens[i] << std::endl;
    }

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
    this_parent_ref = "";
    bool headers = false;
    bool readgroup_present = false;
    while(!headers) {
        std::getline(ifs, line);

        if(line.empty()) {
            return;
        }

        if(line[0] == '@') {
            if(line.substr(0, 3) == "@SQ") {
                std::stringstream ss_sq;
                std::string sq_part;
                ss_sq.str(line);
                std::getline(ss_sq, sq_part, '\t');
                std::getline(ss_sq, sq_part, '\t');
                reference_name = sq_part.substr(3);
                std::getline(ss_sq, sq_part);
                std::string this_len_part = sq_part.substr(3);
                if(!_args.db_names_file.empty()) {
                    if(!_args.rev_db_parent_map.count(reference_name)) {
                        std::cerr << "ERROR: <reference_db>.names file present, but this reference was not detected ";
                        std::cerr << "in the <reference_db>.names file, provided: " << reference_name << std::endl;
                        std::exit(EXIT_FAILURE);
                    }
                    if(!this_parent_ref.empty()) {
                        if(_args.rev_db_parent_map.at(reference_name) != this_parent_ref) {
                            std::cerr << "ERROR: Multiple reference contigs detected that belong to different parent";
                            std::cerr << " relationships. Reads must be aligned to contigs belonging to either a ";
                            std::cerr << "single reference genome or a genome with multiple contigs/segments, whose ";
                            std::cerr << "relations are defined in <reference_db>.names file (see documentation).";
                            std::cerr << "SAM file: " << sam_filepath << std::endl;
                            std::exit(EXIT_FAILURE);
                        }
                    }
                    else {
                        this_parent_ref = _args.rev_db_parent_map.at(reference_name);
                    }
                }
                else {
                    if(!this_parent_ref.empty()) {
                        std::cerr << "ERROR: Multiple reference genomes are only supported if the <reference_db>.names";
                        std::cerr << " file is also present, which defines relations between parent organisms and ";
                        std::cerr << "their children chromosomes/segments (see documentation). SAM file: ";
                        std::cerr << sam_filepath << std::endl;
                        std::exit(EXIT_FAILURE);
                    }
                    this_parent_ref = reference_name;
                }
                this_children_ref.push_back(reference_name);
                ref_lens.push_back(std::stol(this_len_part.c_str()));
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

    for(int i = 0; i < this_children_ref.size(); ++i) {
        if(ref_lens[i] <= 0) {
            std::cerr << "ERROR: Reference sequence length is not positive, provided: ";
            std::cerr << this_children_ref[i] << ", " << ref_lens[i] << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }


    for(int i = 0; i < this_children_ref.size(); ++i) {
        nucleotide_counts[this_children_ref[i]] = std::vector< std::vector< long > >(_iupac_map.size(),
                                                                                     std::vector< long >(ref_lens[i], 0));
        qual_sums[this_children_ref[i]] = std::vector< std::vector< long > >(_iupac_map.size(),
                                                                             std::vector< long >(ref_lens[i], 0));
        mapq_sums[this_children_ref[i]] = std::vector< std::vector< long > >(_iupac_map.size(),
                                                                             std::vector< long >(ref_lens[i], 0));
    }

    std::vector< std::string > res;
    int sam_flag;
    //      0          1           2            3     4    5     6
    // < sam flag, ref name, start pos 1-idx, mapq, cigar, seq, qual >
    res = _parseSamLine(line);
    if((res.size() == 0) || (res[0].empty())) {
        return;
    }
    sam_flag = std::stoi(res[0].c_str());
    if(((sam_flag & 4) == 0) and ((sam_flag & 256) == 0) and ((sam_flag & 2048) == 0)) {
        // Primary alignment
        _addAlignedRead(res[1], res[4], res[5], res[6], std::stol(res[2].c_str()), std::stoi(res[3].c_str()));
    }

    while(std::getline(ifs, line)) {
        res = _parseSamLine(line);
        sam_flag = std::stoi(res[0].c_str());
        if(((sam_flag & 4) == 0) and ((sam_flag & 256) == 0) and ((sam_flag & 2048) == 0)) {
            // Primary alignment
            _addAlignedRead(res[1], res[4], res[5], res[6], std::stol(res[2].c_str()), std::stoi(res[3].c_str()));
        }
    }

    _writePositionalData();

//    printInfo();

    for(auto &[ref, nucl] : nucleotide_counts) {
        while(!_buffer_q->tryPush(sam_sampleid, ref, nucl, qual_sums.at(ref), mapq_sums.at(ref))) {}
    }
}


void ParserJob::_addAlignedRead(const std::string &ref,
                                const std::string &cigar,
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
                    if(!_iupac_map.count(seq[read_idx])) {
                        read_idx++;
                        target_idx++;
                        continue;
                    }
                    nucleotide_counts.at(ref)[_iupac_map.at(seq[read_idx])][target_idx]++;
                    qual_sums.at(ref)[_iupac_map.at(seq[read_idx])][target_idx] += int(qual[read_idx]) - 33;  // Phred 33
                    mapq_sums.at(ref)[_iupac_map.at(seq[read_idx])][target_idx] += mapq;
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
    //      0          1           2            3     4    5     6
    // < sam flag, ref name, start pos 1-idx, mapq, cigar, seq, qual >
    std::vector< std::string > ret;
    std::stringstream this_ss;
    this_ss.str(sam_line);
    std::string this_entry;
    std::getline(this_ss, this_entry, '\t');  // 0. read name
    std::getline(this_ss, this_entry, '\t');  // 1. sam flag
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');  // 2. ref name
    ret.push_back(this_entry);
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

    ofs << "Reference:Index\tA_count,C_count,G_count,T_count\tA_avg_qual,C_avg_qual,G_avg_qual,T_avg_qual\t";
    ofs << "A_avg_mapq,C_avg_mapq,G_avg_mapq,T_avg_mapq" << std::endl;

    for(int r = 0; r < this_children_ref.size(); ++r) {
        std::string ref = this_children_ref[r];
        for(int j = 0; j < ref_lens[r]; ++j) {
            ofs << ref << ':' << (j + 1) << "\t" << nucleotide_counts.at(ref)[0][j];
            for(int i = 1; i < _iupac_map.size(); ++i) {
                ofs << "," << nucleotide_counts.at(ref)[i][j];
            }

            if(nucleotide_counts.at(ref)[0][j] > 0) {
                ofs << "\t" << ((double)qual_sums.at(ref)[0][j] / (double)nucleotide_counts.at(ref)[0][j]);
            }
            else {
                ofs << "\t0";
            }

            for(int i = 1; i < _iupac_map.size(); ++i) {
                if(nucleotide_counts.at(ref)[i][j] > 0) {
                    ofs << "," << ((double)qual_sums.at(ref)[i][j] / (double)nucleotide_counts.at(ref)[i][j]);
                }
                else {
                    ofs << ",0";
                }
            }

            if(nucleotide_counts.at(ref)[0][j] > 0) {
                ofs << "\t" << ((double)mapq_sums.at(ref)[0][j] / (double)nucleotide_counts.at(ref)[0][j]);
            }
            else {
                ofs << "\t0";
            }

            for(int i = 1; i < _iupac_map.size(); ++i) {
                if(nucleotide_counts.at(ref)[i][j] > 0) {
                    ofs << "," << ((double)mapq_sums.at(ref)[i][j] / (double)nucleotide_counts.at(ref)[i][j]);
                }
                else {
                    ofs << ",0";
                }
            }
            ofs << std::endl;
        }
    }

    ofs.close();
}
