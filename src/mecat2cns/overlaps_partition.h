#ifndef OVERLAPS_PARTITION_H
#define OVERLAPS_PARTITION_H

#include <string>	// string
#include <utility>	// pair<>
#include <vector>	// vector<>

#include "../common/alignment.h"

void generate_partition_index_file_name(const std::string& input_file_name, std::string& ret);

void generate_partition_file_name(const std::string& input_file_name, const idx_t part, std::string& ret);

void partition_m4records(const char* m4_file_name, double min_cov_ratio, size_t batch_size, int min_read_size, int num_files);

void partition_candidates(const std::string& input, const std::string& pac_prefix, size_t batch_size, int num_files, idx_t num_reads);

void partition_candidates_reorder(const std::string& input, size_t batch_size, int num_files, idx_t num_reads, const std::vector<idx_t>& read_order, const std::vector<int>& align_counts);

void load_partition_files_info(const char* idx_file_name, std::vector<std::string>& file_info_vec);

void make_read_sort_order(const std::string& input, const std::string& old_pac_prefix, const std::string& sort_file_name, const std::string& pac_prefix, idx_t num_reads, int min_cov, std::vector<idx_t>& read_order, std::vector<std::pair<idx_t, idx_t> >& read_index, std::vector<int>& align_counts);

#endif // OVERLAPS_PARTITION_H
