#include "concurrent_buffer_queue.h"
#include <sstream>
#include <fstream>
#include <algorithm>


ConcurrentBufferQueue::ConcurrentBufferQueue(const int &max_elements) : _max_size(max_elements)
{

}


ConcurrentBufferQueue::~ConcurrentBufferQueue()
{

}


void ConcurrentBufferQueue::run()
{
    std::string output_line, data_line, barcode;
    std::stringstream ss;
    while(!all_jobs_enqueued) {
        while(!tryPop(output_line) && !all_jobs_consumed) {}
        if(!all_jobs_consumed) {

        }
    }
    while(tryPop(output_line)) {

    }

    work_completed = true;


}


bool ConcurrentBufferQueue::tryPush(const std::vector< std::string > &lines,
                                    const std::string &barcode,
                                    const long &reads_processed,
                                    const long &reads_aligned)
{
    std::unique_lock< std::mutex > lock(_mtx);
    if(_q.size() > _max_size) {
        return false;
    }
    for(int i = 0; i < lines.size(); ++i) {
        _q.push(lines[i]);
    }

    aligned_reads_processed[barcode] += reads_aligned;
    total_reads_processed[barcode] += reads_processed;
    return true;
}


bool ConcurrentBufferQueue::tryPop(std::string &item)
{
    std::unique_lock< std::mutex > lock(_mtx);
    if(_q.empty()) {
        if(all_jobs_enqueued) {
            all_jobs_consumed = true;
        }
        return false;
    }

    item = _q.front();
    _q.pop();

    return true;
}


bool ConcurrentBufferQueue::pushHeader(const std::string &barcode, const std::string &header)
{
    std::unique_lock< std::mutex > lock(_mtx);
    if(_headers.count(barcode)) {
        return true;
    }
    _headers[barcode] = header;
    return true;
}
