#include "parser_job.h"
#include <fstream>
#include <sstream>
#include <cassert>


ParserJob::ParserJob(const std::string &parameter_string, ConcurrentBufferQueue* buffer_q) : _buffer_q(buffer_q)
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
                    " reference contig was found in SAM file: " << sam_filepath << std::endl;
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
    printInfo();

//    std::vector< std::string > res;
//    int sam_flag;
//    res = _parseSamLine(line);
//    if((res.size() == 0) || (res[0].empty())) {
//        return;
//    }
//    if(!seen_headers.count(res[0])) {
//        reads_processed++;
//        seen_headers.insert(res[0]);
//    }
//    sam_flag = std::stoi(res[1].c_str());
//    if((sam_flag & 4) == 0) {
//        if(_select) {
//            if(res[2] == genome_select) {
//                contents.push_back(barcode + '|' + line);
//                if(!aligned_headers.count(res[0])) {
//                    reads_aligned++;
//                    aligned_headers.insert(res[0]);
//                }
//            }
//        }
//        else {
//            contents.push_back(barcode + '|' + line);
//            if(!aligned_headers.count(res[0])) {
//                reads_aligned++;
//                aligned_headers.insert(res[0]);
//            }
//        }
//    }
//
//    while(std::getline(ifs, line)) {
//        res = _parseSamLine(line);
//        if(!seen_headers.count(res[0])) {
//            reads_processed++;
//            seen_headers.insert(res[0]);
//        }
//        sam_flag = std::stoi(res[1].c_str());
//        if((sam_flag & 4) == 0) {
//            int temp = sam_flag & 4;
//            if(_select) {
//                if(res[2] == genome_select) {
//                    contents.push_back(barcode + '|' + line);
//                    if(!aligned_headers.count(res[0])) {
//                        reads_aligned++;
//                        aligned_headers.insert(res[0]);
//                    }
//                }
//            }
//            else {
//                contents.push_back(barcode + '|' + line);
//                if(!aligned_headers.count(res[0])) {
//                    reads_aligned++;
//                    aligned_headers.insert(res[0]);
//                }
//            }
//        }
//    }
//
//    while(!_buffer_q->tryPush(contents, barcode, reads_processed, reads_aligned)) {}
}


std::vector< std::string > ParserJob::_parseSamLine(const std::string &sam_line)
{
    std::vector< std::string > ret;
    std::stringstream this_ss;
    this_ss.str(sam_line);
    std::string this_entry;
    std::getline(this_ss, this_entry, '\t');
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');
    ret.push_back(this_entry);
    std::getline(this_ss, this_entry, '\t');
    ret.push_back(this_entry);
    return ret;
}
