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
    while((!all_jobs_enqueued) && (!all_jobs_consumed)) {}

    work_completed = true;
}


bool ConcurrentBufferQueue::tryPush(const std::string &samplename,
                                    const std::vector< std::vector< int > > &nucleotide_counts,
                                    const std::vector< std::vector< long > > &qual_sums,
                                    const std::vector< std::vector < long > > &mapq_sums)
{
    std::unique_lock< std::mutex > lock(_mtx);
    if(all_nucleotide_counts.at(samplename)) {
        std::cerr << "ERROR: Duplicate sample name detected, " << samplename << std::endl;
        std::exit(EXIT_FAILURE);
    }

    all_nucleotide_counts[samplename] = nucleotide_counts;
    all_qual_sums[samplename] = qual_sums;
    all_mapq_sums[samplename] = mapq_sums;

    return true;
}
