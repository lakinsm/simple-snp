#include "concurrent_buffer_queue.h"
#include <sstream>
#include <fstream>
#include <algorithm>


ConcurrentBufferQueue::ConcurrentBufferQueue()
{

}


ConcurrentBufferQueue::~ConcurrentBufferQueue()
{

}


void ConcurrentBufferQueue::run()
{
    std::string output_line, data_line, barcode;
    std::stringstream ss;
    while((!all_jobs_enqueued) || (!all_jobs_consumed)) {}

    work_completed = true;
}


bool ConcurrentBufferQueue::tryPush(const std::string &sample_name,
                                    const std::string &ref_name,
                                    const std::vector< std::vector< int > > &nucleotide_counts,
                                    const std::vector< std::vector< long > > &qual_sums,
                                    const std::vector< std::vector < long > > &mapq_sums)
{
    std::unique_lock< std::mutex > lock(_mtx);
    if(!all_nucleotide_counts.count(sample_name)) {
        all_nucleotide_counts[sample_name];
        all_qual_sums[sample_name];
        all_mapq_sums[sample_name];
    }
    if(!all_nucleotide_counts.at(sample_name).count(ref_name)) {
        std::cerr << "ERROR: Duplicate sample + reference combination detected: " << sample_name;
        std::cerr << ", " << ref_name << std::endl;
        std::exit(EXIT_FAILURE);
    }

    all_nucleotide_counts.at(sample_name)[ref_name] = nucleotide_counts;
    all_qual_sums.at(sample_name)[ref_name] = qual_sums;
    all_mapq_sums.at(sample_name)[ref_name] = mapq_sums;

    return true;
}
