#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <queue>
#include <utility>
#include "args.h"
#include "dispatch_queue.h"
#include "concurrent_buffer_queue.h"
#include "file_finder.h"
#include "fasta_parser.h"


int main(int argc, const char *argv[]) {
    Args args(argc, argv);

    // Load SAM file paths
    FileFinder file_finder;
    std::vector< std::string > sam_files = file_finder.findSamFiles(args.sam_file_dir);

    // Load FASTA reference genome
    FastaParser fasta_parser(args.reference_path);
    fasta_parser.parseFasta();

    DispatchQueue* output_buffer_dispatcher = new DispatchQueue(1, false);
    DispatchQueue* job_dispatcher = new DispatchQueue(args.threads - 1, true);
    ConcurrentBufferQueue* concurrent_q = new ConcurrentBufferQueue();
    output_buffer_dispatcher->dispatch([concurrent_q] () {concurrent_q->run();});

    for(int i = 0; i < sam_files.size(); ++i) {
        std::string this_sam_fp = sam_files[i];
        std::size_t pos1 = this_sam_fp.find_last_of('/');
        std::string this_filename = this_sam_fp.substr(pos1 + 1);
        std::size_t pos2 = this_filename.find_first_of('.');
        std::string this_samplename = this_filename.substr(0, pos2);
        std::string this_param_string = this_sam_fp + '|' + this_samplename;

        while(concurrent_q->num_active_jobs > (args.threads - 2)) {}

        std::unique_ptr< ParserJob > job = std::make_unique< ParserJob > (this_param_string, args.output_dir, concurrent_q);
        job_dispatcher->dispatch(std::move(job));
        concurrent_q->num_active_jobs += 1;
    }

    while(concurrent_q->num_completed_jobs != sam_files.size()) {}
    concurrent_q->all_jobs_enqueued = true;

    while(!concurrent_q->work_completed) {}

    // Each worker thread has written a file with positional counts and info for each sample.  This section is for
    // variant calling across all samples using the thresholds/options specified in args.
    std::ofstream ofs(args.output_dir + "/all_sample_variants.tsv");
    std::vector< std::string > ordered_sample_names;
    for(auto &x : concurrent_q->all_nucleotide_counts) {
        ordered_sample_names.push_back(x.first);
    }
    std::sort(ordered_sample_names.begin(), ordered_sample_names.end());

    ofs << "Position";
    for(int i = 0; i < ordered_sample_names.size(); ++i) {
        ofs << '\t' << ordered_sample_names[i];
    }
    ofs << std::endl;

    std::string this_nucleotides = "ACGT";
    std::map< long, std::map< std::string, std::string > > variant_calls;
    for(int j = 0; j < fasta_parser.seq.size(); ++j) {
        long population_depth = 0;
        // <A, C, G, T>
        std::vector< long > population_allele_counts(4, 0);

        // First pass to look at population metrics
        for(int i = 0; i < population_allele_counts.size(); ++i) {
            for(auto &x : concurrent_q->all_nucleotide_counts) {
                population_depth += x.second[i][j];
                population_allele_counts[i] += x.second[i][j];
            }
        }

        if(population_depth < args.min_inter_sample_depth) {
            continue;
        }

        bool meets_population_threshold = false;
        for(int i = 0; i < population_allele_counts.size(); ++i) {
            meets_population_threshold |= (population_allele_counts[i] > args.min_inter_sample_alt);
        }

        if(!meets_population_threshold) {
            continue;
        }

        // Second pass to establish variants present and their codes
        bool position_has_variant = false;
        std::string alts_present_at_pos = "";
        for(auto &x : concurrent_q->all_nucleotide_counts) {
            long sample_depth = 0;
            for(int i = 0; i < population_allele_counts.size(); ++i) {
                sample_depth += x.second[i][j];
            }

            for(int i = 0; i < population_allele_counts.size(); ++i) {
                double this_allele_freq = (double)x.second[i][j] / (double)sample_depth;
                if((this_allele_freq >= args.min_minor_freq) && (x.second[i][j] >= args.min_intra_sample_alt)) {
                    if(this_nucleotides[i] != fasta_parser.seq[j]) {
                        alts_present_at_pos += this_nucleotides[i];
                        position_has_variant = true;
                    }
                }
            }
        }

        if(!position_has_variant) {
            continue;
        }

        // Third pass to assign variants
        std::map< std::string, std::string > positional_variants;
        for(auto &x : concurrent_q->all_nucleotide_counts) {
            long sample_depth = 0;
            for(int i = 0; i < population_allele_counts.size(); ++i) {
                sample_depth += x.second[i][j];
            }

            if(sample_depth < args.min_intra_sample_depth) {
                std::string low_depth_info = "./.:" + std::to_string(sample_depth) + ":.,.:.,.:.,.";
                positional_variants.insert({x.first, low_depth_info});
                continue;
            }

            std::priority_queue< std::pair< double, std::string > > q;
            for(int i = 0; i < population_allele_counts.size(); ++i) {
                double this_allele_freq = (double)x.second[i][j] / (double)sample_depth;
                if((this_allele_freq >= args.min_minor_freq) && (x.second[i][j] >= args.min_intra_sample_alt)) {
                    std::string var_info;
                    if(this_nucleotides[i] == fasta_parser.seq[j]) {
                        // Reference allele
                        var_info = "0,";
                    }
                    else {
                        std::size_t found = alts_present_at_pos.find(this_nucleotides[i]);
                        var_info = std::to_string(found + 1);
                        var_info += ",";
                    }

                    var_info += std::to_string(x.second[i][j]);
                    var_info += ",";
                    var_info += std::to_string((double)concurrent_q->all_qual_sums.at(x.first)[i][j] / (double)x.second[i][j]);
                    var_info += ",";
                    var_info += std::to_string((double)concurrent_q->all_mapq_sums.at(x.first)[i][j] / (double)x.second[i][j]);
                    q.emplace(this_allele_freq, var_info);
                }
            }
            if(q.size() > 2) {
                std::cerr << "Tri-allelic site detected at sample:position, " << x.first << ":" << (j+1) << std::endl;
                while(!q.empty()) {
                    std::pair< double, std::string > temp_var_info = q.top();
                    std::cout << temp_var_info.first << '\t' << temp_var_info.second << std::endl;
                    q.pop();
                }
                std::cout << std::endl;
                std::exit(EXIT_FAILURE);
            }

            std::string final_var_info = "";
            if(q.size() == 2) {
                std::pair< double, std::string > top_var_info1 = q.top();
                q.pop();
                std::pair< double, std::string > top_var_info2 = q.top();
                std::stringstream ss1, ss2;

                ss1.str(top_var_info1.second);
                ss2.str(top_var_info2.second);

                std::string temp1, temp2;
                std::getline(ss1, temp1, ',');
                std::getline(ss2, temp2, ',');
                final_var_info += temp1 + "/" + temp2 + ":";

                final_var_info += std::to_string(sample_depth) + ":";

                std::getline(ss1, temp1, ',');
                std::getline(ss2, temp2, ',');
                final_var_info += temp1 + "," + temp2 + ":";

                std::getline(ss1, temp1, ',');
                std::getline(ss2, temp2, ',');
                final_var_info += temp1 + "," + temp2 + ":";

                std::getline(ss1, temp1, ',');
                std::getline(ss2, temp2, ',');
                final_var_info += temp1 + "," + temp2;
            }
            else if(q.size() == 1) {
                std::pair< double, std::string > top_var_info = q.top();
                std::stringstream ss;
                ss.str(top_var_info.second);
                std::string temp;
                std::getline(ss, temp, ',');
                final_var_info += temp + "/" + temp + ":";
                final_var_info += std::to_string(sample_depth) + ":";
                std::getline(ss, temp, ',');
                final_var_info += temp + "," + temp + ":";
                std::getline(ss, temp, ',');
                final_var_info += temp + "," + temp + ":";
                std::getline(ss, temp, ',');
                final_var_info += temp + "," + temp;
            }
            else {
                std::string low_depth_info = "./.:" + std::to_string(sample_depth) + ":.,.:.,.:.,.";
                positional_variants.insert({x.first, low_depth_info});
                continue;
            }
            positional_variants.insert({x.first, final_var_info});
        }
        ofs << (j + 1);
        for(int i = 0; i < ordered_sample_names.size(); ++i) {
            ofs << '\t' << positional_variants.at(ordered_sample_names[i]);
        }
        ofs << std::endl;
    }

    ofs.close();
    delete job_dispatcher;
    delete concurrent_q;
    delete output_buffer_dispatcher;

    return 0;
}
