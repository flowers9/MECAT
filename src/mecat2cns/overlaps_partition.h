#ifndef OVERLAPS_PARTITION_H
#define OVERLAPS_PARTITION_H

#include <string>	// string
#include <utility>	// pair<>
#include <vector>	// vector<>

#include "../common/alignment.h"

void generate_partition_index_file_name(const char* m4_file_name, std::string& ret);

void generate_partition_file_name(const char* m4_file_name, const idx_t part, std::string& ret);

void partition_m4records(const char* m4_file_name, double min_cov_ratio, idx_t batch_size, int min_read_size, int num_files);

void partition_candidates(const char* input, idx_t batch_size, int min_read_size, int num_files, idx_t num_reads = 0);

void partition_candidates_reorder(const std::string& input, idx_t batch_size, int num_files, const std::vector<idx_t>& read_order);

void load_partition_files_info(const char* idx_file_name, std::vector<std::string>& file_info_vec);

void make_read_sort_order(const std::string& input, const std::string& sort_file_name, const std::string& pac_prefix, idx_t num_reads, int min_size, int min_cov, std::vector<idx_t>& read_order, std::vector<std::pair<idx_t, idx_t> >& read_info);


#endif // OVERLAPS_PARTITION_H
