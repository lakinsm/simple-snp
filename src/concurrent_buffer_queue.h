#ifndef ASFFAST_CONCURRENT_BUFFER_QUEUE_H
#define ASFFAST_CONCURRENT_BUFFER_QUEUE_H

#include <cassert>
#include <iostream>
#include <mutex>
#include <queue>
#include <atomic>
#include <string>
#include <unordered_map>
#include <vector>


class ConcurrentBufferQueue {
public:
    ConcurrentBufferQueue();
    ~ConcurrentBufferQueue();

    void run();
    bool tryPush(const std::string &sample_name,
                 const std::string &ref_name,
                 const std::vector< std::vector< int > > &nucleotide_counts,
                 const std::vector< std::vector< long > > &qual_sums,
                 const std::vector< std::vector < long > > &mapq_sums,
                 const std::unordered_map< long, std::unordered_map< int, std::vector< long > > > &insertions,
                 const std::unordered_map< long, std::unordered_map< int, std::vector< long > > > &deletions);

    std::atomic< bool > all_jobs_enqueued = ATOMIC_VAR_INIT(false);
    std::atomic< bool > all_jobs_consumed = ATOMIC_VAR_INIT(false);
    std::atomic< bool > work_completed = ATOMIC_VAR_INIT(false);
    std::atomic< int > num_active_jobs = ATOMIC_VAR_INIT(0);
    std::atomic< int > num_completed_jobs = ATOMIC_VAR_INIT(0);

    std::unordered_map< std::string, std::unordered_map< std::string, std::vector< std::vector< int > > > > all_nucleotide_counts;
    std::unordered_map< std::string, std::unordered_map< std::string, std::vector< std::vector< long > > > > all_qual_sums;
    std::unordered_map< std::string, std::unordered_map< std::string, std::vector< std::vector< long > > > > all_mapq_sums;
    std::unordered_map< std::string, std::unordered_map< std::string, std::unordered_map< long, std::unordered_map< int, std::vector< long > > > > > all_insertions;
    std::unordered_map< std::string, std::unordered_map< std::string, std::unordered_map< long, std::unordered_map< int, std::vector< long > > > > > all_deletions;

private:
    std::mutex _mtx;
};


#endif //ASFFAST_CONCURRENT_BUFFER_QUEUE_H
