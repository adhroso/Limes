//
//  Serial.h
//  Limes
//
//  Created by Andi Dhroso on 10/30/13.
//  Copyright (c) 2013 Andi Dhroso. All rights reserved.
//

#ifndef __Limes__Serial__
#define __Limes__Serial__

#include <fstream>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <assert.h>
#include "Types.h"
#include "Timing.hpp"

typedef std::vector<std::pair<std::string::size_type, int> > Lime;
extern Lime limeObjs_target;
extern Lime limeObjs_query;
extern Vector2D vec2D;
extern scottgs::Timing timer;

void remove_invalid_lime_candidates(Candidates &candidates, const RawData &d1, const RawData &d2);
void remove_invalid_lime_candidates_reverse(Candidates &candidates, const RawData &query_data, const RawData &target_data);
void find_candidate_limes(Chromosome &query_chr, Candidates &candidates);
void find_limes_in_forward_direction(Chromosome &query_chr, const RawData &query_data, const RawData &target_data, Candidates &candidates);
void find_limes_in_reverse_direction(Chromosome &query_chr, const RawData &query_data, const RawData &target_data, Candidates &candidates);
void find_all_limes_with_pair(const RawData &target_data, const RawData &query_data, const Chromosome& target_chr);

void find_all_limes_serial_approach(const std::string &path_to_genome_a, const std::string &path_to_genome_b, const std::string &ext, const std::string &output, const std::string &progressFile);
void find_all_limes_serial_approach(const std::string &path_to_genome_a, const std::string &path_to_genome_b, const std::string &ext, const std::string &output);
void find_all_limes_with_sequence(const std::string &path_to_sequences_a, const std::string &path_to_genome_b, const std::string &ext, const std::string &output);
void generate_limes(const std::vector<std::string> &g1, const std::vector<std::string> &g2, const std::string &limes, const std::string &progress);
void generate_limes(const std::string &seqFile, const std::vector<std::string> &g, const std::string &lime, const std::string &progress);
void generate_limes(const std::string &fileG1, const std::string &fileG2, const std::string &limes, const std::string &progress_file);

void find_all_limes_between_unassembled_genomes_cluster(const std::string &path_to_file, const std::string &path_to_genome, const std::string &ext, const std::string &output, const std::string &progress);
void find_all_limes_between_unassembled_genomes_non_cluster(const std::string &path_to_genome_a, const std::string &path_to_genome_b, const std::string &ext, const std::string &output, const std::string &progress);


void reverse_complement(RawData &);

bool file_exists(const std::string &path_to_file);
void load_next_batch(std::ifstream &in, std::vector<std::pair<std::string, std::string> > &sequences, const int batch_size=1000);

#endif /* defined(__Limes__Serial__) */
