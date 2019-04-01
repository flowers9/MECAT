#ifndef OPTIONS_H
#define OPTIONS_H

#include "../common/defs.h"
#include <string>

#define INPUT_TYPE_CAN 	0
#define INPUT_TYPE_M4	1

struct ConsensusOptions
{
    int         input_type;
    const char* m4;
    const char* reads;
    const char* corrected_reads;
    const char* grid_options;
    const char* grid_options_split;
    int         num_threads;
    idx_t       batch_size;
    double      min_mapping_ratio;
    int         min_align_size;
    int         min_cov;
    idx_t       min_size;
    bool        print_usage_info;
    int         tech;
    int         num_partition_files;
    int         job_index;
    int         reads_to_correct;
    int		grid_start_delay;
    int		full_reads;
    idx_t	read_buffer_size;
};

void
print_usage(const char* prog);

int
parse_arguments(int argc, char* argv[], ConsensusOptions& t);

void
print_options(ConsensusOptions& t);

std::string
make_options(const ConsensusOptions& t);

typedef ConsensusOptions ReadsCorrectionOptions;

#endif // OPTIONS_H
