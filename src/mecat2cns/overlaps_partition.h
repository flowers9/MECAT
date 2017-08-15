#ifndef OVERLAPS_PARTITION_H
#define OVERLAPS_PARTITION_H

#include <vector>

#include "../common/alignment.h"

void 
generate_partition_index_file_name(const char* m4_file_name, std::string& ret);

void
generate_partition_file_name(const char* m4_file_name, const idx_t part, std::string& ret);

void
partition_m4records(const char* m4_file_name, const double min_cov_ratio, const idx_t batch_size, const int min_read_size);

void
partition_candidates(const char* input, const idx_t batch_size, const int min_read_size);

struct PartitionFileInfo
{
    std::string file_name;
    idx_t min_seq_id;
    idx_t max_seq_id;
};

void
load_partition_files_info(const char* idx_file_name, std::vector<PartitionFileInfo>& file_info_vec);

#endif // OVERLAPS_PARTITION_H
