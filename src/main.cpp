#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <queue>
#include <utility>
#include <cmath>
#include "args.h"
#include "dispatch_queue.h"
#include "concurrent_buffer_queue.h"
#include "file_finder.h"
#include "fasta_parser.h"
#include "vcf_writer.h"


int main(int argc, const char *argv[]) {
    Args args(argc, argv);

    // Optionally load database annotations
    if(args.db_ann_file != "") {
        std::ifstream ifs6(args.db_ann_file, std::ios::in);
        std::string ann_line, ann_acc, ann_entry;
        std::stringstream ann_ss;
        std::getline(ifs6, ann_line);  // Skip header
        while(std::getline(ifs6, ann_line)) {
            if(ann_line.empty()) {
                continue;
            }
            ann_ss.clear();
            ann_ss.str(ann_line);
            std::getline(ann_ss, ann_acc, ',');
            if(!args.db_ann_map.count(ann_acc)) {
                args.db_ann_map[ann_acc] = std::vector< std::vector< std::string > >(1, std::vector< std::string >());
            }
            else {
                args.db_ann_map.at(ann_acc).push_back(std::vector< std::string >());
            }
            int j_vec_idx = args.db_ann_map.at(ann_acc).size() - 1;
            for(int i = 0; i < 5; ++i) {
                std::getline(ann_ss, ann_entry, ',');
                args.db_ann_map.at(ann_acc)[j_vec_idx].push_back(ann_entry);
            }
        }
        ifs6.close();
    }

    // Optionally load database names
    if(!args.db_names_file.empty()) {
        std::ifstream ifs7(args.db_names_file, std::ios::in);
        std::string names_line, names_parent, names_child, names_alias;
        std::stringstream names_ss;
        while(std::getline(ifs7, names_line)) {
            if(names_line.empty()) {
                continue;
            }
            names_ss.clear();
            names_ss.str(names_line);
            std::getline(names_ss, names_parent, ',');
            std::getline(names_ss, names_alias, ',');
            std::size_t div_pos = names_parent.find(':');
            if(div_pos == std::string::npos) {
                // No children are present
                if(args.db_parent_name_map.count(names_parent)) {
                    std::cerr << "ERROR: Parent chromosomes must be unique if no children are present,";
                    std::cerr << " (duplicate detected): ";
                    std::cerr << names_parent << std::endl;
                    exit(EXIT_FAILURE);
                }
                args.db_parent_name_map[names_parent] = names_alias;
                args.db_child_name_map[names_parent] = names_alias;
                args.db_parent_map[names_parent] = std::vector< std::string > {names_parent};
                args.rev_db_parent_map[names_parent] = names_parent;
            }
            else {
                // Children are present
                // names_alias format: parent_name\tchild_name
                names_child = names_parent.substr(div_pos + 1);
                names_parent.erase(div_pos);
                if(args.db_parent_map.count(names_parent)) {
                    args.db_parent_map.at(names_parent).push_back(names_child);
                }
                else {
                    args.db_parent_map[names_parent] = std::vector< std::string > {names_child};
                    if(args.rev_db_parent_map.count(names_child)) {
                        std::cerr << "ERROR: Child chromosomes must be unique (duplicate detected): ";
                        std::cerr << names_parent << " -> " << names_child << std::endl;
                        exit(EXIT_FAILURE);
                    }
                    args.rev_db_parent_map[names_child] = names_parent;
                }
                if(args.db_child_name_map.count(names_child)) {
                    std::cerr << "ERROR: Child chromosomes must be unique (duplicate detected): ";
                    std::cerr << names_parent << " -> " << names_child << ": " << names_alias << std::endl;
                    exit(EXIT_FAILURE);
                }
                if(!args.db_parent_name_map.count(names_parent)) {
                    std::size_t parent_found = names_alias.find_last_of(',');
                    args.db_parent_name_map[names_parent] = names_alias.substr(0, parent_found - 1);
                }
                std::size_t child_found = names_alias.find_last_of(',');
                args.db_child_name_map[names_child] = names_alias.substr(child_found + 1);
            }
        }
        ifs7.close();
    }

    // Load SAM file paths
    FileFinder file_finder;
    std::vector< std::string > sam_files = file_finder.findSamFiles(args.sam_file_dir);

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

        std::unique_ptr< ParserJob > job = std::make_unique< ParserJob > (this_param_string, args.output_dir, concurrent_q, args);
        job_dispatcher->dispatch(std::move(job));
        concurrent_q->num_active_jobs += 1;
    }
    concurrent_q->all_jobs_enqueued = true;

    while(concurrent_q->num_completed_jobs != sam_files.size()) {}
    concurrent_q->all_jobs_consumed = true;

    while(!concurrent_q->work_completed) {}

    // Each worker thread has written a file with positional counts and info for each sample.  This section is for
    // variant calling across all samples using the thresholds/options specified in args.

    // Check to ensure all SAM files have a valid reference (parent/child relationship for multi-chromosome refs)
    // Store ordered vectors of samples/refs
    std::vector< std::string > ordered_sample_names;
    std::string this_parent_ref = "";
    std::string sample_ref;
    for(auto &[sample, ref_map] : concurrent_q->all_nucleotide_counts) {
        for(auto &[ref, nucl] : ref_map) {
            if(!args.db_names_file.empty()) {
                sample_ref = args.rev_db_parent_map.at(ref);
            }
            else {
                sample_ref = ref;
            }
            if(this_parent_ref.empty()) {
                this_parent_ref = sample_ref;
            }
            else {
                if(this_parent_ref != sample_ref) {
                    std::cerr << "ERROR: All SAM files must be aligned to the same reference. If the reference used ";
                    std::cerr << "has multiple chromosomes/segments, they must be defined in <reference_db>.names ";
                    std::cerr << "(see documentation). Clashing references: " << this_parent_ref << ", " << sample_ref;
                    std::cerr << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            }
        }
        ordered_sample_names.push_back(sample);
    }
    std::sort(ordered_sample_names.begin(), ordered_sample_names.end());

    std::vector< std::string > ordered_refs;
    for(int i = 0; i < args.db_parent_map.at(this_parent_ref).size(); ++i) {
        ordered_refs.push_back(args.db_parent_map.at(this_parent_ref)[i]);
    }

    std::ofstream ofs(args.output_dir + "/all_sample_variants.tsv");
    std::ofstream ofs2(args.output_dir + "/dominant_population_variants.tsv");

    // VCF Writer
    std::string command_string = "simple_snp " + args.sam_file_dir + " " + args.output_dir + " " + args.reference_path;
    command_string += " -t " + std::to_string(args.threads) + " -a " + std::to_string(args.min_intra_sample_alt);
    command_string += " -A " + std::to_string(args.min_inter_sample_alt) + " -d " + std::to_string(args.min_intra_sample_depth);
    command_string += " -D " + std::to_string(args.min_inter_sample_depth) + " -f " + std::to_string(args.min_minor_freq);
    command_string += " -F " + std::to_string(args.min_major_freq);
    if(!args.db_ann_file.empty()) {
        command_string += " -n " + args.db_ann_file;
    }

    // Load FASTA reference genome
    FastaParser fasta_parser(args.reference_path);
    fasta_parser.parseFasta(ordered_refs);

    std::vector< long > contig_lens;
    for(int i = 0; i < ordered_refs.size(); ++i) {
        contig_lens.push_back(fasta_parser.headers_lens.at(ordered_refs[i]));
    }

    std::string vcf_path = args.output_dir + "/dominant_population_variants.vcf";
    VcfWriter vcf_writer(vcf_path);
    vcf_writer.open();
    vcf_writer.writeHeaders(args.reference_path,
                            command_string,
                            ordered_refs,
                            contig_lens);
    vcf_writer.writeSamples(ordered_sample_names);

    ofs << "Position";
    ofs2 << "Position";
    for(int i = 0; i < ordered_sample_names.size(); ++i) {
        ofs << '\t' << ordered_sample_names[i];
        ofs2 << '\t' << ordered_sample_names[i];
    }
    ofs << std::endl;
    ofs2 << std::endl;

    std::string this_nucleotides = "ACGT";
    for(int r = 0; r < ordered_refs.size(); ++r) {
        std::string this_ref = ordered_refs[r];
        std::string this_seq = fasta_parser.headers_seqs.at(this_ref);
        for(long j = 0; j < this_seq.length(); ++j) {
            long population_depth = 0;
            // <A, C, G, T>
            std::vector< long > population_allele_counts(4, 0);
            std::unordered_map< int, long > population_insertions;
            std::unordered_map< int, long > population_deletions;

            // First pass to look at population metrics
            for(auto &[sample, ref_map] : concurrent_q->all_nucleotide_counts) {
                std::vector< std::vector< int > > *nucl = &ref_map.at(this_ref);
                std::unordered_map< long, std::unordered_map< int, std::vector< long > > > *ins = &concurrent_q->all_insertions.at(sample).at(this_ref);
                std::unordered_map< long, std::unordered_map< int, std::vector< long > > > *del = &concurrent_q->all_deletions.at(sample).at(this_ref);
                long sample_depth = 0;
                for(int i = 0; i < population_allele_counts.size(); ++i) {
                    population_depth += (*nucl)[i][j];
                    sample_depth += (*nucl)[i][j];
                }

                if((*ins).count(j)) {
                    for(auto &[len, ins_vec] : (*ins).at(j)) {
                        population_depth += ins_vec[0];
                        sample_depth += ins_vec[0];
                    }
                }

                if((*del).count(j)) {
                    for(auto &[len, del_vec] : (*del).at(j)) {
                        population_depth += del_vec[0];
                        sample_depth += del_vec[0];
                    }
                }

                if(population_depth < args.min_inter_sample_depth) {
                    continue;
                }

                for(int i = 0; i < population_allele_counts.size(); ++i) {
                    double this_allele_freq = (double)(*nucl)[i][j] / (double)sample_depth;
                    if(this_allele_freq >= args.min_major_freq) {
                        if(this_nucleotides[i] != this_seq[j]) {
                            population_allele_counts[i] += (*nucl)[i][j];
                        }
                    }
                }

                // indel frequency is calculated across all indel lengths to identify candidate indels at a given
                // position. This is to mitigate the effect of nanopore sequencing noise, particular when indels
                // are identified near homopolymer runs.
                if((*ins).count(j)) {
                    double this_ins_freq = 0;
                    for(auto &[len, ins_vec] : (*ins).at(j)) {
                         this_ins_freq += (double)ins_vec[0];
                    }
                    this_ins_freq /= (double)sample_depth;
                    if(this_ins_freq >= args.min_major_freq) {
                        for(auto &[len, ins_vec] : (*ins).at(j)) {
                            if(!population_insertions.count(len)) {
                                population_insertions[len] = ins_vec[0];
                            }
                            else {
                                population_insertions.at(len) += ins_vec[0];
                            }
                        }
                    }
                }

                if((*del).count(j)) {
                    double this_del_freq = 0;
                    for(auto &[len, del_vec] : (*del).at(j)) {
                        this_del_freq += (double)del_vec[0];
                    }
                    this_del_freq /= (double)sample_depth;
                    if(this_del_freq >= args.min_major_freq) {
                        for(auto &[len, del_vec] : (*del).at(j)) {
                            if(!population_deletions.count(len)) {
                                population_deletions[len] = del_vec[0];
                            }
                            else {
                                population_deletions[len] += del_vec[0];
                            }
                        }
                    }
                }
            }

            bool meets_population_threshold = false;
            for(int i = 0; i < population_allele_counts.size(); ++i) {
                meets_population_threshold |= (population_allele_counts[i] > args.min_inter_sample_alt);
            }

            long population_ins_sums = 0;
            for(auto &[len, val] : population_insertions) {
                population_ins_sums += val;
            }
            meets_population_threshold |= (population_ins_sums > args.min_inter_sample_alt);

            long population_del_sums = 0;
            for(auto &[len, val] : population_deletions) {
                population_del_sums += val;
            }
            meets_population_threshold |= (population_del_sums > args.min_inter_sample_alt);

            if(!meets_population_threshold) {
                continue;
            }

            // Second pass to establish variants present and their codes
            vcfLineData vcf_line_data;
            vcf_line_data.dp = 0;

            bool position_has_variant = false;
            bool position_has_major_variant = false;
            std::string alts_present_at_pos = "";
            for(auto &[sample, ref_map] : concurrent_q->all_nucleotide_counts) {
                std::vector< std::vector< int > > *nucl = &ref_map.at(this_ref);
                std::unordered_map< long, std::unordered_map< int, std::vector< long > > > *ins = &concurrent_q->all_insertions.at(sample).at(this_ref);
                std::unordered_map< long, std::unordered_map< int, std::vector< long > > > *del = &concurrent_q->all_deletions.at(sample).at(this_ref);
                long sample_depth = 0;
                for(int i = 0; i < population_allele_counts.size(); ++i) {
                    sample_depth += (*nucl)[i][j];
                }

                if((*ins).count(j)) {
                    for(auto &[len, ins_vec] : (*ins).at(j)) {
                        sample_depth += ins_vec[0];
                    }
                }

                if((*del).count(j)) {
                    for(auto &[len, del_vec] : (*del).at(j)) {
                        sample_depth += del_vec[0];
                    }
                }

                vcf_line_data.dp += sample_depth;

                for(int i = 0; i < population_allele_counts.size(); ++i) {
                    double this_allele_freq = (double)(*nucl)[i][j] / (double)sample_depth;
                    if((this_allele_freq >= args.min_minor_freq) && ((*nucl)[i][j] >= args.min_intra_sample_alt) && (sample_depth > args.min_intra_sample_depth)) {
                        if(this_nucleotides[i] != this_seq[j]) {
                            if(alts_present_at_pos.find(this_nucleotides[i]) == std::string::npos) {
                                alts_present_at_pos += this_nucleotides[i];
                            }
                            position_has_variant = true;
                            if(this_allele_freq >= args.min_major_freq) {
                                position_has_major_variant = true;
                            }
                        }
                    }
                }

                if((*ins).count(j)) {
                    double this_ins_freq = 0;
                    long this_ins_count = 0;
                    for(auto &[len, ins_vec] : (*ins).at(j)) {
                        this_ins_count += ins_vec[0];
                        this_ins_freq += (double)ins_vec[0];
                    }
                    this_ins_freq /= (double)sample_depth;
                    if((this_ins_freq >= args.min_minor_freq) && (this_ins_count >= args.min_intra_sample_alt) && (sample_depth > args.min_intra_sample_depth)) {
                        std::cout << this_ref << ':' << std::to_string(j+1) << "\tInsertion\t" << this_ins_count;
                        std::cout << '\t' << this_ins_freq << std::endl;
                        for(auto &[len, ins_vec] : (*ins).at(j)) {
                            std::cout << '\t' << len << '\t' << ins_vec[0] << '\t' << ins_vec[1] << '\t';
                            std::cout << ins_vec[2] << '\t' << ins_vec[3] << std::endl;
                        }
                        if(alts_present_at_pos.find('I') == std::string::npos) {
                            alts_present_at_pos += 'I';
                        }
                        position_has_variant = true;
                        if(this_ins_freq >= args.min_major_freq) {
                            position_has_major_variant = true;
                        }
                    }
                }

                if((*del).count(j)) {
                    double this_del_freq = 0;
                    long this_del_count = 0;
                    for(auto &[len, del_vec] : (*del).at(j)) {
                        this_del_count += del_vec[0];
                        this_del_freq += (double)del_vec[0];
                    }
                    this_del_freq /= (double)sample_depth;
                    if((this_del_freq >= args.min_minor_freq) && (this_del_count >= args.min_intra_sample_alt) && (sample_depth > args.min_intra_sample_depth)) {
                        std::cout << this_ref << ':' << std::to_string(j+1) << "\tDeletion\t" << this_del_count;
                        std::cout << '\t' << this_del_freq << std::endl;
                        for(auto &[len, del_vec] : (*del).at(j)) {
                            std::cout << '\t' << len << '\t' << del_vec[0] << '\t' << del_vec[1] << '\t';
                            std::cout << del_vec[2] << std::endl;
                        }
                        if(alts_present_at_pos.find('D') == std::string::npos) {
                            alts_present_at_pos += 'D';
                        }
                        position_has_variant = true;
                        if(this_del_freq >= args.min_major_freq) {
                            position_has_major_variant = true;
                        }
                    }
                }
            }

            if(!position_has_variant) {
                continue;
            }

            vcf_line_data.chrom = this_ref;
            vcf_line_data.ref = this_seq[j];
            vcf_line_data.pos = j+1;
            vcf_line_data.qual = 0;
            vcf_line_data.ns = 0;
            vcf_line_data.ro = 0;
            vcf_line_data.mqmr = 0;
            vcf_line_data.ao_sum = 0;
            vcf_line_data.nsa = 0;
            for(int i = 0; i < alts_present_at_pos.size(); ++i) {

                if(alts_present_at_pos[i] == 'I') {
                    vcf_line_data.type.push_back("ins");
                }
                else if(alts_present_at_pos[i] == 'D') {
                    vcf_line_data.type.push_back("del");
                }
                else {
                    vcf_line_data.type.push_back("snp");
                }

                // TODO:  needs to be moved below with incorporation of indels
                vcf_line_data.cigar.push_back("1X");
                vcf_line_data.af.push_back(0);
                vcf_line_data.alt.push_back("");
                vcf_line_data.alt[i] += alts_present_at_pos[i];
                vcf_line_data.ao.push_back(0);
                vcf_line_data.mqm.push_back(0);
                vcf_line_data.alt_ns.push_back(0);
                vcf_line_data.ac.push_back(0);
            }

            // Third pass to assign variants
            std::map< std::string, std::string > positional_variants;
            std::map< std::string, std::string > vcf_variants;
            for(auto &[sample, ref_map] : concurrent_q->all_nucleotide_counts) {
                std::vector< std::vector< int > > *nucl = &ref_map.at(this_ref);
                std::vector< std::vector< long > > *qual = &concurrent_q->all_qual_sums.at(sample).at(this_ref);
                std::vector< std::vector< long > > *mapq = &concurrent_q->all_mapq_sums.at(sample).at(this_ref);
                std::unordered_map< long, std::unordered_map< int, std::vector< long > > > *ins = &concurrent_q->all_insertions.at(sample).at(this_ref);
                std::unordered_map< long, std::unordered_map< int, std::vector< long > > > *del = &concurrent_q->all_deletions.at(sample).at(this_ref);
                long sample_depth = 0;
                int ref_allele_count;
                double ref_qual;
                for(int i = 0; i < population_allele_counts.size(); ++i) {
                    sample_depth += (*nucl)[i][j];
                    if(this_nucleotides[i] == this_seq[j]) {
                        ref_allele_count = (*nucl)[i][j];
                        ref_qual = (double)(*qual)[i][j];
                        vcf_line_data.mqmr += (double)(*mapq)[i][j];
                    }
                }

                if((*ins).count(j)) {
                    for(auto &[len, ins_vec] : (*ins).at(j)) {
                        sample_depth += ins_vec[0];
                    }
                }

                if((*del).count(j)) {
                    for(auto &[len, del_vec] : (*del).at(j)) {
                        sample_depth += del_vec[0];
                    }
                }

                if(sample_depth < args.min_intra_sample_depth) {
                    std::string low_depth_info = "./.:" + std::to_string(sample_depth) + ":.:.:.";
                    positional_variants.insert({sample, low_depth_info});
                    // GT:DP:AD:RO:QR:AO:QA
                    std::string low_vcf_info = "./.:" + std::to_string(sample_depth) + ":.";
                    for(int i = 0; i < alts_present_at_pos.size(); ++i) {
                        low_vcf_info += ",.";
                    }
                    low_vcf_info += ":.:.";
                    for(int i = 1; i < alts_present_at_pos.size(); ++i) {
                        low_vcf_info += ",.";
                    }
                    low_vcf_info += ":.:.";
                    for(int i = 1; i < alts_present_at_pos.size(); ++i) {
                        low_vcf_info += ",.";
                    }
                    vcf_variants.insert({sample, low_vcf_info});
                    continue;
                }

//            std::cout << '\t' << sample << '\t' << j+1 << std::flush;

                std::priority_queue< std::pair< double, std::string > > q;
                for(int i = 0; i < population_allele_counts.size(); ++i) {
                    double this_allele_freq = (double)(*nucl)[i][j] / (double)sample_depth;
                    if((this_allele_freq >= args.min_minor_freq) && ((*nucl)[i][j] >= args.min_intra_sample_alt)) {
                        std::string var_info;
                        if(this_nucleotides[i] == this_seq[j]) {
                            // Reference allele
                            var_info = "0,";
                        }
                        else {
                            std::size_t found = alts_present_at_pos.find(this_nucleotides[i]);
                            var_info = std::to_string(found + 1);
                            var_info += ",";
                            vcf_line_data.mqm[found] += (double)(*mapq)[i][j];
                            vcf_line_data.ao[found] += (*nucl)[i][j];
                            vcf_line_data.ao_sum += (*nucl)[i][j];
                            vcf_line_data.qual += (double)(*qual)[i][j];
                        }

                        var_info += std::to_string((*nucl)[i][j]);
                        var_info += ",";
                        var_info += std::to_string((double)(*qual)[i][j] / (double)(*nucl)[i][j]);
                        var_info += ",";
                        var_info += std::to_string((double)(*mapq)[i][j] / (double)(*nucl)[i][j]);
                        var_info += ",";
                        var_info += std::to_string(ref_allele_count);
                        var_info += ",";
                        if(ref_allele_count > 0) {
                            var_info += std::to_string(ref_qual / (double)ref_allele_count);
                        }
                        else {
                            var_info += ".";
                        }
                        q.emplace(this_allele_freq, var_info);
                        vcf_line_data.ro += ref_allele_count;
                    }
                }

                if((*ins).count(j)) {
                    double this_ins_freq = 0;
                    long this_ins_count = 0;
                    for(auto &[len, ins_vec] : (*ins).at(j)) {
                        this_ins_count += ins_vec[0];
                        this_ins_freq += (double)ins_vec[0];
                    }
                    this_ins_freq /= (double)sample_depth;
                    if((this_ins_freq >= args.min_minor_freq) && (this_ins_count >= args.min_intra_sample_alt)) {
//                        std::cout << this_ref << ':' << std::to_string(j+1) << "\tInsertion\t" << this_ins_count;
//                        std::cout << '\t' << this_ins_freq << std::endl;
//                        for(auto &[len, ins_vec] : (*ins).at(j)) {
//                            std::cout << '\t' << len << '\t' << ins_vec[0] << '\t' << ins_vec[1] << '\t';
//                            std::cout << ins_vec[2] << '\t' << ins_vec[3] << std::endl;
//                        }
                    }
                }

                if((*del).count(j)) {
                    double this_del_freq = 0;
                    long this_del_count = 0;
                    for(auto &[len, del_vec] : (*del).at(j)) {
                        this_del_count += del_vec[0];
                        this_del_freq += (double)del_vec[0];
                    }
                    this_del_freq /= (double)sample_depth;
                    if((this_del_freq >= args.min_minor_freq) && (this_del_count >= args.min_intra_sample_alt)) {
//                        std::cout << this_ref << ':' << std::to_string(j+1) << "\tDeletion\t" << this_del_count;
//                        std::cout << '\t' << this_del_freq << std::endl;
//                        for(auto &[len, del_vec] : (*del).at(j)) {
//                            std::cout << '\t' << len << '\t' << del_vec[0] << '\t' << del_vec[1] << '\t';
//                            std::cout << del_vec[2] << std::endl;
//                        }
                    }
                }


//                if(q.size() > 2) {
//                    std::cerr << "Tri-allelic site detected at sample:position, " << sample << " ";
//                    std::cerr << this_ref << ":" << (j+1) << std::endl;
//                    while(!q.empty()) {
//                        std::pair< double, std::string > temp_var_info = q.top();
//                        std::cout << temp_var_info.first << '\t' << temp_var_info.second << std::endl;
//                        q.pop();
//                    }
//                    std::cout << std::endl;
//                    std::exit(EXIT_FAILURE);
//                }
//
//                std::string final_var_info = "";
//                std::string final_vcf_info = "";
//                if(q.size() == 2) {
//                    std::pair< double, std::string > top_var_info1 = q.top();
//                    q.pop();
//                    std::pair< double, std::string > top_var_info2 = q.top();
//                    std::stringstream ss1, ss2;
//
//                    ss1.str(top_var_info1.second);
//                    ss2.str(top_var_info2.second);
//
//                    std::string temp1, temp2;
//                    std::string ro, qr;
//                    std::string gt1, ao1, gq1, qa1;
//                    std::string gt2, ao2, gq2, qa2;
//
//                    // Genotype
//                    std::getline(ss1, gt1, ',');
//                    std::getline(ss2, gt2, ',');
//                    final_var_info += gt1 + "/" + gt2 + ":";
//                    final_vcf_info += gt1 + "/" + gt2 + ":";
//
//                    if((gt1 != "0") or (gt2 != "0")) {
//                        vcf_line_data.nsa++;
//                    }
//                    vcf_line_data.ns++;
//
//                    // Depth
//                    final_var_info += std::to_string(sample_depth) + ":";
//                    final_vcf_info += std::to_string(sample_depth) + ":";
//
//                    // Allele count
//                    std::getline(ss1, ao1, ',');
//                    std::getline(ss2, ao2, ',');
//                    final_var_info += ao1 + "," + ao2 + ":";
//
//                    // Mean quality score
//                    std::getline(ss1, qa1, ',');
//                    std::getline(ss2, qa2, ',');
//                    final_var_info += temp1 + "," + temp2 + ":";
//
//                    // Mean mapq score
//                    std::getline(ss1, temp1, ',');
//                    std::getline(ss2, temp2, ',');
//                    final_var_info += temp1 + "," + temp2 + ":";
//
//                    if(gt1 != "0") {
//                        vcf_line_data.alt_ns[std::stoi(gt1.c_str()) - 1] += 1;
//                    }
//                    if(gt2 != "0") {
//                        vcf_line_data.alt_ns[std::stoi(gt2.c_str()) - 1] += 1;
//                    }
//
//                    // Ref allele count
//                    std::getline(ss1, ro, ',');
//                    final_var_info += ro + ":";
//
//                    // Mean ref allele qual score
//                    std::getline(ss1, qr, ',');
//                    final_var_info += qr;
//
//                    std::vector< int > sample_vcf_ao(vcf_line_data.ao.size(), 0);
//                    std::vector< double > sample_vcf_qa(vcf_line_data.ao.size(), 0);
//
//                    int gt1_idx = std::stoi(gt1.c_str()) - 1;
//                    int gt2_idx = std::stoi(gt2.c_str()) - 1;
//                    if(gt1_idx >= 0) {
//                        sample_vcf_ao[gt1_idx] = std::stoi(ao1.c_str());
//                        sample_vcf_qa[gt1_idx] = std::stod(qa1.c_str());
//                        vcf_line_data.ac[gt1_idx]++;
//                    }
//                    if(gt2_idx >= 0) {
//                        sample_vcf_ao[gt2_idx] = std::stoi(ao2.c_str());
//                        sample_vcf_qa[gt2_idx] = std::stod(qa2.c_str());
//                        vcf_line_data.ac[gt2_idx]++;
//                    }
//
//                    final_vcf_info += ro;
//                    for(int i = 0; i < sample_vcf_ao.size(); ++i) {
//                        final_vcf_info += ',' + std::to_string(sample_vcf_ao[i]);
//                    }
//
//                    final_vcf_info +=  ":" + ro + ":" + qr + ":";
//                    final_vcf_info += std::to_string(sample_vcf_ao[0]);
//                    for(int i = 1; i < sample_vcf_ao.size(); ++i) {
//                        final_vcf_info += ',' + std::to_string(sample_vcf_ao[i]);
//                    }
//                    final_vcf_info += ":" + std::to_string(sample_vcf_qa[0]);
//                    for(int i = 1; i < sample_vcf_qa.size(); ++i) {
//                        final_vcf_info += ',' + std::to_string(sample_vcf_qa[i]);
//                    }
//                }
//                else if(q.size() == 1) {
//                    std::pair< double, std::string > top_var_info = q.top();
//                    std::stringstream ss;
//                    ss.str(top_var_info.second);
//                    std::string temp;
//                    std::string gt, ro, ao, gq, qr, qa;
//                    std::getline(ss, gt, ',');
//                    final_var_info += gt + "/" + gt + ":";
//                    final_vcf_info += gt + "/" + gt + ":";
//                    final_var_info += std::to_string(sample_depth) + ":";
//                    final_vcf_info += std::to_string(sample_depth) + ":";
//                    std::getline(ss, ao, ',');
//                    final_var_info += ao + "," + ao + ":";
//                    std::getline(ss, qa, ',');
//                    final_var_info += qa + "," + qa + ":";
//                    std::getline(ss, temp, ',');
//                    final_var_info += temp + "," + temp + ":";
//                    std::getline(ss, ro, ',');
//                    final_var_info += ro + ":";
//                    std::getline(ss, qr, ',');
//                    final_var_info += qr;
//
//                    if(gt == "0") {
//                        int sample_nucl_idx;
//                        final_vcf_info += ro;
//                        for(int i = 0; i < vcf_line_data.ao.size(); ++i) {
//                            sample_nucl_idx = this_nucleotides.find(alts_present_at_pos[i]);
//                            final_vcf_info += ',' + std::to_string((*nucl)[sample_nucl_idx][j]);
//                        }
//                        final_vcf_info += ":" + ro + ":" + qr + ":";
//                        sample_nucl_idx = this_nucleotides.find(alts_present_at_pos[0]);
//                        final_vcf_info += std::to_string((*nucl)[sample_nucl_idx][j]);
//                        for(int i = 1; i < vcf_line_data.ao.size(); ++i) {
//                            sample_nucl_idx = this_nucleotides.find(alts_present_at_pos[i]);
//                            final_vcf_info += ',' + std::to_string((*nucl)[sample_nucl_idx][j]);
//                        }
//                        sample_nucl_idx = this_nucleotides.find(alts_present_at_pos[0]);
//                        if((*nucl)[sample_nucl_idx][j] > 0) {
//                            final_vcf_info += ":" + std::to_string((double)(*qual)[sample_nucl_idx][j] /
//                                                                   (double)(*nucl)[sample_nucl_idx][j]);
//                        }
//                        else {
//                            final_vcf_info += ":.";
//                        }
//
//                        for(int i = 1; i < vcf_line_data.ao.size(); ++i) {
//                            sample_nucl_idx = this_nucleotides.find(alts_present_at_pos[i]);
//                            final_vcf_info += ',';
//                            if((*nucl)[sample_nucl_idx][j] > 0) {
//                                final_vcf_info += std::to_string((double)(*qual)[sample_nucl_idx][j] /
//                                                                 (double)(*nucl)[sample_nucl_idx][j]);
//                            }
//                            else {
//                                final_vcf_info += '.';
//                            }
//                        }
//                    }
//                    else {
//                        final_vcf_info += ro + ',' + ao + ":" + ro + ":" + qr + ":" + ao + ":" + qa;
//                    }
//
//                    vcf_line_data.ns++;
//                    if(gt != "0") {
//                        int gt_idx = std::stoi(gt.c_str()) - 1;
//                        vcf_line_data.ac[gt_idx] += 2;
//                        vcf_line_data.nsa++;
//                    }
//                }
//                else {
//                    std::string low_depth_info = "./.:" + std::to_string(sample_depth) + ":.:.:.";
//                    positional_variants.insert({sample, low_depth_info});
//                    std::string low_vcf_info = "./.:" + std::to_string(sample_depth) + ":.";
//                    for(int i = 0; i < alts_present_at_pos.size(); ++i) {
//                        low_vcf_info += ",.";
//                    }
//                    low_vcf_info += ":.:.";
//                    for(int i = 1; i < alts_present_at_pos.size(); ++i) {
//                        low_vcf_info += ",.";
//                    }
//                    low_vcf_info += ":.:.";
//                    for(int i = 1; i < alts_present_at_pos.size(); ++i) {
//                        low_vcf_info += ",.";
//                    }
//                    vcf_variants.insert({sample, low_vcf_info});
//                    continue;
//                }
//                positional_variants.insert({sample, final_var_info});
//                vcf_variants.insert({sample, final_vcf_info});
            }
//            ofs << this_ref << ':' << (j + 1);
//            for(int i = 0; i < ordered_sample_names.size(); ++i) {
//                ofs << '\t' << positional_variants.at(ordered_sample_names[i]);
//            }
//            ofs << std::endl;
//
//            if(position_has_major_variant) {
//                ofs2 << this_ref << ':' << (j + 1);
//                for(int i = 0; i < ordered_sample_names.size(); ++i) {
//                    ofs2 << '\t' << positional_variants.at(ordered_sample_names[i]);
//                }
//                ofs2 << std::endl;
//            }
//
//            vcf_line_data.qual = std::log((double)vcf_line_data.ao_sum) * (vcf_line_data.qual / (double)vcf_line_data.ao_sum);
//            for(int i = 0; i < vcf_line_data.alt.size(); ++i) {
//                vcf_line_data.af[i] = (double)vcf_line_data.ao[i] / (double)vcf_line_data.dp;
//                vcf_line_data.mqm[i] /= (double)vcf_line_data.ao[i];
//            }
//            vcf_line_data.mqmr /= (double)vcf_line_data.ro;
//
//            vcf_writer.writeSampleData(vcf_line_data, vcf_variants);
        }
    }

    ofs.close();
    ofs2.close();
    vcf_writer.close();
    delete job_dispatcher;
    delete concurrent_q;
    delete output_buffer_dispatcher;

    return 0;
}
