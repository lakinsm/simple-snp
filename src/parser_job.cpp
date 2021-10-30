#include "parser_job.h"
#include <fstream>
#include <sstream>


ParserJob::ParserJob(const std::string &parameter_string, ConcurrentBufferQueue* buffer_q) : _buffer_q(buffer_q)
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
    std::cout << smaplename << '\t' << sam_filepath << std::endl;
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
    while(!headers) {
        std::getline(ifs, line);

        if(line.empty()) {
            return;
        }

        if(line[0] == '@') {
            if(line.substr(0, 3) == "@RG") {
                std::sstream ss_rg;
                std::string rg_part;
                ss_rg.str(line);
                std::getline(ss_rg, rg_part, '\t');
                std::getline(ss_rg, rg_part, '\t');
                sam_readgroup = rg_part.substr(3);
                std::getline(ss_rg, rg_part);
                sam_sampleid = rg_part.substr(3);
            }
        }
        else {
            headers = true;
        }
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
