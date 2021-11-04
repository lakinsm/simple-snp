#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
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
    std::ofstream ofs2(args.output_dir + "/dominant_population_variants.tsv");

    std::vector< std::string > ordered_sample_names;
    for(auto &x : concurrent_q->all_nucleotide_counts) {
        ordered_sample_names.push_back(x.first);
    }
    std::sort(ordered_sample_names.begin(), ordered_sample_names.end());

    // VCF Writer
    std::string command_string = "simple_snp " + args.sam_file_dir + " " + args.output_dir + " " + args.reference_path;
    command_string += " -t " + std::to_string(args.threads) + " -a " + std::to_string(args.min_intra_sample_alt);
    command_string += " -A " + std::to_string(args.min_inter_sample_alt) + " -d " + std::to_string(args.min_intra_sample_depth);
    command_string += " -D " + std::to_string(args.min_inter_sample_depth) + " -f " + std::to_string(args.min_minor_freq);
    command_string += " -F " + std::to_string(args.min_major_freq);
    std::vector< std::string > contig_names = {fasta_parser.header};
    std::vector< long > contig_lens = {(long)fasta_parser.seq.size()};
    std::string vcf_path = args.output_dir + "/dominant_population_variants.vcf";
    VcfWriter vcf_writer(vcf_path);
    vcf_writer.open();
    vcf_writer.writeHeaders(args.reference_path,
                            command_string,
                            contig_names,
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
    std::map< long, std::map< std::string, std::string > > variant_calls;
    for(int j = 0; j < fasta_parser.seq.size(); ++j) {
        long population_depth = 0;
        // <A, C, G, T>
        std::vector< long > population_allele_counts(4, 0);

        // First pass to look at population metrics
        for(auto &x : concurrent_q->all_nucleotide_counts) {
            long sample_depth = 0;
            for(int i = 0; i < population_allele_counts.size(); ++i) {
                population_depth += x.second[i][j];
                sample_depth += x.second[i][j];
            }

            for(int i = 0; i < population_allele_counts.size(); ++i) {
                double this_allele_freq = (double)x.second[i][j] / (double)sample_depth;
                if(this_allele_freq >= args.min_major_freq) {
                    if(this_nucleotides[i] != fasta_parser.seq[j]) {
                        population_allele_counts[i] += x.second[i][j];
                    }
                }
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
        vcfLineData vcf_line_data;
        vcf_line_data.dp = 0;

        bool position_has_variant = false;
        bool position_has_major_variant = false;
        std::string alts_present_at_pos = "";
        for(auto &x : concurrent_q->all_nucleotide_counts) {
            long sample_depth = 0;
            for(int i = 0; i < population_allele_counts.size(); ++i) {
                sample_depth += x.second[i][j];
            }

            vcf_line_data.dp += sample_depth;

            for(int i = 0; i < population_allele_counts.size(); ++i) {
                double this_allele_freq = (double)x.second[i][j] / (double)sample_depth;
                if((this_allele_freq >= args.min_minor_freq) && (x.second[i][j] >= args.min_intra_sample_alt) && (sample_depth > args.min_intra_sample_depth)) {
                    if(this_nucleotides[i] != fasta_parser.seq[j]) {
                        if(alts_present_at_pos.find(this_nucleotides[i]) == std::string::npos) {
                            alts_present_at_pos += this_nucleotides[i];
                        }
                        position_has_variant = true;
                    }
                }
                if((this_allele_freq >= args.min_major_freq) && (x.second[i][j] >= args.min_intra_sample_alt) && (sample_depth > args.min_intra_sample_depth)) {
                    if(this_nucleotides[i] != fasta_parser.seq[j]) {
                        std::cout << x.first << '\t' << (j+1) << '\t' << '\t' << this_nucleotides[i];
                        std::cout << '\t' << fasta_parser.seq[j] << '\t';
                        std::cout << this_allele_freq << '\t' << x.second[i][j] << std::endl;
                        position_has_major_variant = true;
                    }
                }
            }
        }

        if(!position_has_variant) {
            continue;
        }

        vcf_line_data.chrom = fasta_parser.header;
        vcf_line_data.ref = fasta_parser.seq[j];
        vcf_line_data.pos = j+1;
        vcf_line_data.qual = 0;
        vcf_line_data.ns = 0;
        vcf_line_data.ro = 0;
        vcf_line_data.mqmr = 0;
        vcf_line_data.ac = 0;
        vcf_line_data.ao_sum = 0;
        vcf_line_data.nsa = 0;
        for(int i = 0; i < alts_present_at_pos.size(); ++i) {
            vcf_line_data.type.push_back("snp");
            vcf_line_data.cigar.push_back("1X");
            vcf_line_data.af.push_back(0);
            vcf_line_data.alt.push_back("");
            vcf_line_data.alt[i] += alts_present_at_pos[i];
            vcf_line_data.ao.push_back(0);
            vcf_line_data.mqm.push_back(0);
            vcf_line_data.alt_ns.push_back(0);
        }

        // Third pass to assign variants
        std::map< std::string, std::string > positional_variants;
        std::map< std::string, std::string > vcf_variants;
        for(auto &x : concurrent_q->all_nucleotide_counts) {
            long sample_depth = 0;
            int ref_allele_count;
            double ref_qual;
            for(int i = 0; i < population_allele_counts.size(); ++i) {
                sample_depth += x.second[i][j];
                if(this_nucleotides[i] == fasta_parser.seq[j]) {
                    ref_allele_count = x.second[i][j];
                    ref_qual = (double)concurrent_q->all_qual_sums.at(x.first)[i][j];
                    vcf_line_data.mqmr += (double)concurrent_q->all_mapq_sums.at(x.first)[i][j];
                }
            }

            if(sample_depth < args.min_intra_sample_depth) {
                std::string low_depth_info = "./.:" + std::to_string(sample_depth) + ":.:.:.";
                positional_variants.insert({x.first, low_depth_info});
                // GT:DP:AD:RO:QR:AO:QA
                std::string low_vcf_info = "./.:" + std::to_string(sample_depth) + ":.:.:.:.:.";
                vcf_variants.insert({x.first, low_vcf_info});
                continue;
            }

//            std::cout << '\t' << x.first << '\t' << j+1 << std::flush;

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
                        vcf_line_data.mqm[found] += (double)concurrent_q->all_mapq_sums.at(x.first)[i][j];
                        vcf_line_data.ao[found] += x.second[i][j];
                        vcf_line_data.ao_sum += x.second[i][j];
                        vcf_line_data.qual += (double)concurrent_q->all_qual_sums.at(x.first)[i][j];
                    }

                    var_info += std::to_string(x.second[i][j]);
                    var_info += ",";
                    var_info += std::to_string((double)concurrent_q->all_qual_sums.at(x.first)[i][j] / (double)x.second[i][j]);
                    var_info += ",";
                    var_info += std::to_string((double)concurrent_q->all_mapq_sums.at(x.first)[i][j] / (double)x.second[i][j]);
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

//            std::cout << "\tCheck1" << std::flush;

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
            std::string final_vcf_info = "";
            if(q.size() == 2) {
//                std::cout << "\tq2" << std::flush;
                std::pair< double, std::string > top_var_info1 = q.top();
                q.pop();
                std::pair< double, std::string > top_var_info2 = q.top();
                std::stringstream ss1, ss2;

                ss1.str(top_var_info1.second);
                ss2.str(top_var_info2.second);

                std::string temp1, temp2;
                std::string ro, qr;
                std::string gt1, ao1, gq1, qa1;
                std::string gt2, ao2, gq2, qa2;

                // Genotype
                std::getline(ss1, gt1, ',');
                std::getline(ss2, gt2, ',');
                final_var_info += gt1 + "/" + gt2 + ":";
                final_vcf_info += gt1 + "/" + gt2 + ":";

                if((gt1 != "0") or (gt2 != "0")) {
                    vcf_line_data.nsa++;
                }
                vcf_line_data.ns++;

                // Depth
                final_var_info += std::to_string(sample_depth) + ":";
                final_vcf_info += std::to_string(sample_depth) + ":";

                // Allele count
                std::getline(ss1, ao1, ',');
                std::getline(ss2, ao2, ',');
                final_var_info += ao1 + "," + ao2 + ":";

                // Mean quality score
                std::getline(ss1, qa1, ',');
                std::getline(ss2, qa2, ',');
                final_var_info += temp1 + "," + temp2 + ":";

                if(gt1 != "0") {
                    vcf_line_data.ac++;
                }
                if(gt2 != "0") {
                    vcf_line_data.ac++;
                }

                // Mean mapq score
                std::getline(ss1, temp1, ',');
                std::getline(ss2, temp2, ',');
                final_var_info += temp1 + "," + temp2 + ":";

                if(gt1 != "0") {
                    vcf_line_data.alt_ns[std::stoi(gt1.c_str()) - 1] += 1;
                }
                if(gt2 != "0") {
                    vcf_line_data.alt_ns[std::stoi(gt2.c_str()) - 1] += 1;
                }

                // Ref allele count
                std::getline(ss1, ro, ',');
                final_var_info += ro + ":";

                // Mean ref allele qual score
                std::getline(ss1, qr, ',');
                final_var_info += qr;

                std::vector< int > sample_vcf_ao(vcf_line_data.ao.size(), 0);
                std::vector< double > sample_vcf_qa(vcf_line_data.ao.size(), 0);

                int gt1_idx = std::stoi(gt1.c_str()) - 1;
                int gt2_idx = std::stoi(gt2.c_str()) - 1;

//                std::cout << '\t' << vcf_line_data.ao.size() << '\t' << gt1_idx << "/" << gt2_idx << std::endl;

                if(gt1_idx >= 0) {
                    sample_vcf_ao[gt1_idx] = std::stoi(ao1.c_str());
                    sample_vcf_qa[gt1_idx] = std::stod(qa1.c_str());
                }
                if(gt2_idx >= 0) {
                    sample_vcf_ao[gt2_idx] = std::stoi(ao2.c_str());
                    sample_vcf_qa[gt2_idx] = std::stod(qa2.c_str());
                }

                final_vcf_info += ro;
                for(int i = 0; i < sample_vcf_ao.size(); ++i) {
                    final_vcf_info += ',' + std::to_string(sample_vcf_ao[i]);
                }

                final_vcf_info +=  ":" + ro + ":" + qr + ":";
                final_vcf_info += std::to_string(sample_vcf_ao[0]);
                for(int i = 1; i < sample_vcf_ao.size(); ++i) {
                    final_vcf_info += ',' + std::to_string(sample_vcf_ao[i]);
                }
                final_vcf_info += ":" + std::to_string(sample_vcf_qa[0]);
                for(int i = 1; i < sample_vcf_qa.size(); ++i) {
                    final_vcf_info += ',' + std::to_string(sample_vcf_qa[i]);
                }
            }
            else if(q.size() == 1) {
//                std::cout << "\tq1" << std::endl;
                std::pair< double, std::string > top_var_info = q.top();
                std::stringstream ss;
                ss.str(top_var_info.second);
                std::string temp;
                std::string gt, ro, ao, gq, qr, qa;
                std::getline(ss, gt, ',');
                final_var_info += gt + "/" + gt + ":";
                final_vcf_info += gt + "/" + gt + ":";
                final_var_info += std::to_string(sample_depth) + ":";
                final_vcf_info += std::to_string(sample_depth) + ":";
                std::getline(ss, ao, ',');
                final_var_info += ao + "," + ao + ":";
                std::getline(ss, qa, ',');
                final_var_info += qa + "," + qa + ":";
                std::getline(ss, temp, ',');
                final_var_info += temp + "," + temp + ":";
                std::getline(ss, ro, ',');
                final_var_info += ro + ":";
                std::getline(ss, qr, ',');
                final_var_info += qr;

                if(gt == "0") {
                    int sample_nucl_idx;
                    final_vcf_info += ro;
                    for(int i = 0; i < vcf_line_data.ao.size(); ++i) {
                        sample_nucl_idx = this_nucleotides.find(alts_present_at_pos[i]);
                        final_vcf_info += ',' + std::to_string(x.second[sample_nucl_idx][j]);
                    }
                    final_vcf_info += ":" + ro + ":" + qr + ":";
                    sample_nucl_idx = this_nucleotides.find(alts_present_at_pos[0]);
                    final_vcf_info += std::to_string(x.second[sample_nucl_idx][j]);
                    for(int i = 1; i < vcf_line_data.ao.size(); ++i) {
                        sample_nucl_idx = this_nucleotides.find(alts_present_at_pos[i]);
                        final_vcf_info += ',' + std::to_string(x.second[sample_nucl_idx][j]);
                    }
                    sample_nucl_idx = this_nucleotides.find(alts_present_at_pos[0]);
                    if(x.second[sample_nucl_idx][j] > 0) {
                        final_vcf_info += ":" + std::to_string((double)concurrent_q->all_qual_sums.at(x.first)[sample_nucl_idx][j] /
                                                               (double)x.second[sample_nucl_idx][j]);
                    }
                    else {
                        final_vcf_info += ":.";
                    }

                    for(int i = 1; i < vcf_line_data.ao.size(); ++i) {
                        sample_nucl_idx = this_nucleotides.find(alts_present_at_pos[i]);
                        final_vcf_info += ',';
                        if(x.second[sample_nucl_idx][j] > 0) {
                            final_vcf_info += std::to_string((double)concurrent_q->all_qual_sums.at(x.first)[sample_nucl_idx][j] /
                                    (double)x.second[sample_nucl_idx][j]);
                        }
                        else {
                            final_vcf_info += '.';
                        }
                    }
                }
                else {
                    final_vcf_info += ro + ',' + ao + ":" + ro + ":" + qr + ":" + ao + ":" + qa;
                }

                vcf_line_data.ns++;
                if(gt != "0") {
                    vcf_line_data.ac += 2;
                    vcf_line_data.nsa++;
                }
            }
            else {
                std::string low_depth_info = "./.:" + std::to_string(sample_depth) + ":.:.:.";
                positional_variants.insert({x.first, low_depth_info});
                std::string low_vcf_info = "./.:" + std::to_string(sample_depth) + ":.:.:.:.:.";
                vcf_variants.insert({x.first, low_vcf_info});
                continue;
            }
            positional_variants.insert({x.first, final_var_info});
            vcf_variants.insert({x.first, final_vcf_info});
        }
        ofs << (j + 1);
        for(int i = 0; i < ordered_sample_names.size(); ++i) {
            ofs << '\t' << positional_variants.at(ordered_sample_names[i]);
        }
        ofs << std::endl;

        if(position_has_major_variant) {
            ofs2 << (j + 1);
            for(int i = 0; i < ordered_sample_names.size(); ++i) {
                ofs2 << '\t' << positional_variants.at(ordered_sample_names[i]);
            }
            ofs2 << std::endl;
        }

        vcf_line_data.qual = std::log((double)vcf_line_data.ao_sum) * (vcf_line_data.qual / (double)vcf_line_data.ao_sum);
        for(int i = 0; i < vcf_line_data.alt.size(); ++i) {
            vcf_line_data.af[i] = (double)vcf_line_data.ao[i] / (double)vcf_line_data.dp;
            vcf_line_data.mqm[i] /= (double)vcf_line_data.ao[i];
        }
        vcf_line_data.mqmr /= (double)vcf_line_data.ro;

        vcf_writer.writeSampleData(vcf_line_data, vcf_variants);
    }

    ofs.close();
    ofs2.close();
    vcf_writer.close();
    delete job_dispatcher;
    delete concurrent_q;
    delete output_buffer_dispatcher;

    return 0;
}
