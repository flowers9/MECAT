#ifndef OVERLAPS_PARTITION_H
#define OVERLAPS_PARTITION_H

#include <vector>	// vector<>
#include <string>	// string

#include "../common/alignment.h"

void generate_partition_index_file_name(const char* m4_file_name, std::string& ret);

void generate_partition_file_name(const char* m4_file_name, const idx_t part, std::string& ret);

void partition_m4records(const char* m4_file_name, double min_cov_ratio, idx_t batch_size, int min_read_size, int num_files);

void partition_candidates(const char* input, idx_t batch_size, int min_read_size, int num_files, idx_t num_reads = 0);

void load_partition_files_info(const char* idx_file_name, std::vector<std::string>& file_info_vec);

#endif // OVERLAPS_PARTITION_H
