#include "options.h"

#include <cstring>
#include <unistd.h>

#include <iostream>
#include <sstream>
#include <string>

static int input_type_pacbio		= 1;
static int num_threads_pacbio		= 1;
static double mapping_ratio_pacbio	= 0.6;
static int align_size_pacbio		= 1000;
static int cov_pacbio			= 4;
static int min_size_pacbio		= 2000;
static bool print_usage_pacbio		= false;
static int tech_pacbio			= TECH_PACBIO;

static int input_type_nanopore		= 1;
static int num_threads_nanopore		= 1;
static double mapping_ratio_nanopore	= 0.4;
static int align_size_nanopore		= 400;
static int cov_nanopore			= 6;
static int min_size_nanopore		= 2000;
static bool print_usage_nanopore	= false;
static int tech_nanopore		= TECH_NANOPORE;

static int default_tech = TECH_PACBIO;
static int num_partition_files = 0;
static int full_reads = 0;
static idx_t read_buffer_size = 0;
static int preprocess_reads = 0;
static int reorder_reads = 0;
static idx_t batch_size = 8589934592;		// 8 GB

static const char input_type_n    = 'i';
static const char num_threads_n   = 't';
static const char batch_size_n    = 'p';
static const char mapping_ratio_n = 'r';
static const char align_size_n    = 'a';
static const char cov_n           = 'c';
static const char min_size_n      = 'l';
static const char usage_n         = 'h';
static const char tech_n          = 'x';
static const char num_partition_files_n = 'k';
static const char grid_options_n  = 'G';
static const char grid_options_split_n = 'S';
static const char job_index_n     = 'I';
static const char reads_to_correct_n = 'R';
static const char grid_start_delay_n = 'D';
static const char full_reads_n    = 'F';
static const char read_buffer_size_n = 'b';
static const char preprocess_reads_n = 'P';
static const char reorder_reads_n = 'O';

void
print_pacbio_default_options()
{
	std::cerr << "-" << input_type_n << " " << input_type_pacbio
		 << " -" << num_threads_n << " " << num_threads_pacbio
		 << " -" << mapping_ratio_n << " " << mapping_ratio_pacbio
		 << " -" << align_size_n << " " << align_size_pacbio
		 << " -" << cov_n << " " << cov_pacbio << " "
		 << " -" << min_size_n << " " << min_size_pacbio
		 << "\n";
}

void print_nanopore_default_options()
{
	std::cerr << "-" << input_type_n << " " << input_type_nanopore
		 << " -" << num_threads_n << " " << num_threads_nanopore
		 << " -" << mapping_ratio_n << " " << mapping_ratio_nanopore
		 << " -" << align_size_n << " " << align_size_nanopore
		 << " -" << cov_n << " " << cov_nanopore
		 << " -" << min_size_n << " " << min_size_nanopore
		 << "\n";
}

// given options, recreate arguments from the command line
std::string
make_options(const ConsensusOptions& options)
{
        std::ostringstream cmd;
	cmd << " -" << input_type_n << " " << (options.input_type == INPUT_TYPE_CAN ? 0 : 1);
	if (options.num_threads > -1) {
		cmd << " -" << num_threads_n << " " << options.num_threads;
	}
	if (options.batch_size && options.batch_size != batch_size) {
		cmd << " -" << batch_size_n << " " << options.batch_size;
	}
	if (options.min_mapping_ratio >= 0 && ((options.tech == TECH_PACBIO && options.min_mapping_ratio != mapping_ratio_pacbio) || (options.tech != TECH_PACBIO && options.min_mapping_ratio != mapping_ratio_nanopore))) {
		cmd << " -" << mapping_ratio_n << " " << options.min_mapping_ratio;
	}
	if (options.min_align_size >= 0 && ((options.tech == TECH_PACBIO && options.min_align_size != align_size_pacbio) || (options.tech != TECH_PACBIO && options.min_align_size != align_size_nanopore))) {
		cmd << " -" << align_size_n << " " << options.min_align_size;
	}
	if (options.min_cov >= 0 && ((options.tech == TECH_PACBIO && options.min_cov != cov_pacbio) || (options.tech != TECH_PACBIO && options.min_cov != cov_nanopore))) {
		cmd << " -" << cov_n << " " << options.min_cov;
	}
	if (options.min_size >= 0 && ((options.tech == TECH_PACBIO && options.min_size != min_size_pacbio) || (options.tech != TECH_PACBIO && options.min_size != min_size_nanopore))) {
		cmd << " -" << min_size_n << " " << options.min_size;
	}
	if (options.num_partition_files > 0 && options.num_partition_files != num_partition_files) {
		cmd << " -" << num_partition_files_n << " " << options.num_partition_files;
	}
	if (options.grid_options != NULL) {
		cmd << " -" << grid_options_n << " \"" << options.grid_options << "\"";
	}
	if (options.grid_options_split != NULL) {
		cmd << " -" << grid_options_split_n << " \"" << options.grid_options_split << "\"";
	}
	if (options.job_index != -1) {
		cmd << " -" << job_index_n << " " << options.job_index;
	}
	if (options.reads_to_correct) {
		cmd << " -" << reads_to_correct_n << " " << options.reads_to_correct;
	}
	if (options.grid_start_delay) {
		cmd << " -" << grid_start_delay_n << " " << options.grid_start_delay;
	}
	if (options.full_reads) {
		cmd << " -" << full_reads_n;
	}
	if (options.read_buffer_size) {
		cmd << " -" << read_buffer_size_n << " " << options.read_buffer_size;
	}
	if (options.preprocess_reads) {
		cmd << " -" << preprocess_reads_n;
	}
	if (options.reorder_reads) {
		cmd << " -" << reorder_reads_n;
	}
	cmd << " " << options.m4;
	cmd << " " << options.reads;
	cmd << " " << options.corrected_reads;
        return cmd.str();
}

void print_usage(const char* prog) {
	std::cerr << "USAGE:\n"
		<< prog << " [options] input reads output\n"
		<< "\n"
		<< "OPTIONS:\n"
		<< "-" << tech_n << " <0/1>\tsequencing platform: 0 = PACBIO, 1 = NANOPORE [0]\n"
		<< "-" << input_type_n << " <0/1>\tinput type: 0 = candidate, 1 = m4\n"
		<< "-" << num_threads_n << " <Integer>\tnumber of threads (CPUs)\n"
		<< "-" << batch_size_n << " <Integer>\tmemory used for holding candidates [8 GB]\n"
		<< "-" << mapping_ratio_n << " <Real>\tminimum mapping ratio\n"
		<< "-" << align_size_n << " <Integer>\tminimum overlap size\n"
		<< "-" << cov_n << " <Integer>\tminimum coverage under consideration\n"
		<< "-" << min_size_n << " <Integer>\tminimum length of corrected sequence\n"
		<< "-" << num_partition_files_n << " <Integer>\tnumber of partition files when partitioning overlap results (if 0, then use system limit)\n"
		<< "-" << grid_options_n << " <String>\toptions for grid submission\n"
		<< "-" << grid_options_split_n << " <String>\toptions for split grid submission\n"
		<< "-" << reads_to_correct_n << " <Integer>\tnumber of reads to correct [all]\n"
		<< "-" << grid_start_delay_n << " <Integer>\tseconds to delay between starting grid jobs\n"
		<< "-" << full_reads_n << "\t\toutput full reads, not just the corrected parts\n"
		<< "-" << read_buffer_size_n << " <Integer>\tbytes of memory to buffer reads [no buffer]\n"
		<< "-" << preprocess_reads_n << "\t\tconvert reads from fasta to fasta db before processing (implied by -" << read_buffer_size_n << " and -" << reorder_reads_n << ")\n"
		<< "-" << reorder_reads_n << "\t\treorder reads before processing to improve memory efficiency\n"
		<< "-" << usage_n << "\t\tprint usage info\n"
		<< "\n"
		<< "If 'x' is set to be '0' (pacbio), then the other options have the following default values: \n";
	print_pacbio_default_options();
	std::cerr << "\n"
		<< "If 'x' is set to be '1' (nanopore), then the other options have the following default values: \n";
	print_nanopore_default_options();
}

ConsensusOptions init_consensus_options(const int tech) {
	ConsensusOptions t;
	t.m4                    = NULL;
	t.reads                 = NULL;
	t.corrected_reads       = NULL;
	t.grid_options          = NULL;
	t.grid_options_split    = NULL;
	t.num_partition_files	= num_partition_files;
	t.job_index		= -1;
	t.reads_to_correct	= 0;
	t.grid_start_delay	= 0;
	t.full_reads		= full_reads;
	t.read_buffer_size	= read_buffer_size;
	t.preprocess_reads	= preprocess_reads;
	t.reorder_reads		= reorder_reads;
	t.batch_size            = batch_size;
	if (tech == TECH_PACBIO) {
		t.input_type            = input_type_pacbio;
		t.num_threads           = num_threads_pacbio;
		t.min_mapping_ratio     = mapping_ratio_pacbio;
		t.min_align_size        = align_size_pacbio;
		t.min_cov               = cov_pacbio;
		t.min_size              = min_size_pacbio;
		t.print_usage_info      = print_usage_pacbio;
		t.tech                  = tech_pacbio;
	} else {
		t.input_type            = input_type_nanopore;
		t.num_threads           = num_threads_nanopore;
		t.min_mapping_ratio     = mapping_ratio_nanopore;
		t.min_align_size        = align_size_nanopore;
		t.min_cov               = cov_nanopore;
		t.min_size              = min_size_nanopore;
		t.print_usage_info      = print_usage_nanopore;
		t.tech                  = tech_nanopore;
	}
	return t;
}

int detect_tech(int argc, char* argv[]) {
	int t = default_tech;
	char tech_nstr[64]; tech_nstr[0] = '-'; tech_nstr[1] = tech_n; tech_nstr[2] = '\0';
	for (int i = 0; i < argc; ++i) {
		if (strcmp(tech_nstr, argv[i]) == 0) {
			if (i + 1 == argc) {
				fprintf(stderr, "argument to option '%c' is missing.\n", tech_n);
				t = -1;
			}
			std::cout << tech_nstr << "\n";
			if (argv[i + 1][0] == '0') {
				t = TECH_PACBIO;
			} else if (argv[i + 1][0] == '1') {
				t = TECH_NANOPORE;
			} else {
				fprintf(stderr, "invalid argument to option '%c': %s\n", tech_n, argv[i + 1]);
				t = -1;
			}
			break;
		}
	}
	return t;
}

int parse_arguments(int argc, char* argv[], ConsensusOptions& t) {
	bool parse_success = true;
	int tech = detect_tech(argc, argv);
	if (tech == -1) {
		return 1;
	} 
	t = init_consensus_options(tech);
	int opt_char;
	char err_char;
	opterr = 0;
	while ((opt_char = getopt(argc, argv, "i:t:p:r:a:c:l:x:hG:I:n:R:S:D:k:Fb:PO")) != -1) {
		switch (opt_char) {
			case input_type_n:
				if (optarg[0] == '0') {
					t.input_type = INPUT_TYPE_CAN;
				} else if (optarg[0] == '1') {
					t.input_type = INPUT_TYPE_M4;
				} else {
					fprintf(stderr, "invalid argument to option '%c': %s\n", input_type_n, optarg);
					return 1;
				}
				break;
			case num_threads_n:
				t.num_threads = atoi(optarg);
				break;
			case batch_size_n:
				t.batch_size = atoll(optarg);
				break;
			case mapping_ratio_n:
				t.min_mapping_ratio = atof(optarg);
				break;
			case align_size_n:
				t.min_align_size = atoi(optarg);
				break;
			case cov_n:
				t.min_cov = atoi(optarg);
				break;
			case min_size_n:
				t.min_size = atoll(optarg);
				break;
			case usage_n:
				t.print_usage_info = true;
				break;
			case grid_options_n:
				t.grid_options = optarg;
				break;
			case grid_options_split_n:
				t.grid_options_split = optarg;
				break;
			case job_index_n:
				t.job_index = atoi(optarg);
				break;
			case reads_to_correct_n:
				t.reads_to_correct = atoi(optarg);
				break;
			case grid_start_delay_n:
				t.grid_start_delay = atoi(optarg);
				break;
			case tech_n:
				break;
			case num_partition_files_n:
				t.num_partition_files = atoi(optarg);
				break;
			case full_reads_n:
				t.full_reads = 1;
				break;
			case read_buffer_size_n:
				t.read_buffer_size = atoll(optarg);
				// implies preprocess_reads, so fall through here
			case preprocess_reads_n:
				t.preprocess_reads = 1;
				break;
			case reorder_reads_n:
				t.preprocess_reads = 1;
				t.reorder_reads = 1;
				break;
			case '?':
				err_char = (char)optopt;
				fprintf(stderr, "unrecognised option '%c'\n", err_char);
				return 1;
				break;
			case ':':
				err_char = (char)optopt;
				fprintf(stderr, "argument to option '%c' is missing.\n", err_char);
				return 1;
				break;
		}
	}
	if (t.num_threads < 0) {
		std::cerr << "cpu threads must be greater than 0\n";
		parse_success = false;
	}
	if (t.batch_size == 0) {
		std::cerr << "batch size must be greater than 0\n";
		parse_success = false;
	}
	if (t.min_mapping_ratio < 0.0) {
		std::cerr << "mapping ratio must be >= 0.0\n";
		parse_success = false;
	}
	if (t.min_cov < 0) {
		std::cerr << "coverage must be >= 0\n";
		parse_success = false;
	}
	if (t.min_size < 0) {
		std::cerr << "sequence size must be >= 0\n";
		parse_success = false;
	}
	if (t.reads_to_correct < 0) {
		std::cerr << "number of reads must be >= 0\n";
		parse_success = false;
	}
	if (t.grid_start_delay < 0) {
		std::cerr << "grid start delay must be >= 0\n";
		parse_success = false;
	}
	if (t.read_buffer_size < 0) {
		std::cerr << "read buffer size must be greater than 0\n";
		parse_success = false;
	}
	if (argc - optind < 3) {
		return 1;
	}
	t.m4 = argv[argc - 3];
	t.reads = argv[argc - 2];
	t.corrected_reads = argv[argc - 1];
	return parse_success ? 0 : 1;
}

void
print_options(ConsensusOptions& t)
{
	std::cout << "input_type:\t"	<< t.input_type << "\n";
	if (t.m4) std::cout << "reads\t" << t.m4 << "\n";
	if (t.reads) std::cout << "output\t" << t.reads << "\n";
	if (t.corrected_reads) std::cout << "m4\t" << t.corrected_reads << "\n";
	if (t.grid_options) std::cout << "grid\t" << t.grid_options << "\n";
	if (t.grid_options_split) std::cout << "grid_split\t" << t.grid_options_split << "\n";
	if (t.full_reads) std::cout << "full reads\n";
	if (t.read_buffer_size) std::cout << "read_buffer_size\t" << t.read_buffer_size << "\n";
	std::cout << "number of threads:\t" << t.num_threads << "\n";
	std::cout << "batch size:\t" << t.batch_size << "\n";
	std::cout << "mapping ratio:\t" << t.min_mapping_ratio << "\n";
	std::cout << "align size:\t" << t.min_align_size << "\n";
	std::cout << "cov:\t" << t.min_cov << "\n";
	std::cout << "min size:\t" << t.min_size << "\n";
	std::cout << "partition files:\t" << t.num_partition_files << "\n";
	std::cout << "tech:\t" << t.tech << "\n";
}
