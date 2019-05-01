#include "reads_correction_can.h"
#include "reads_correction_m4.h"
#include "overlaps_partition.h"
#include "options.h"
#include "packed_db.h"			// PackedDB

#include <fcntl.h>	// S_IRUSR, S_IXUSR
#include <list>
#include <sstream>
#include <string>	// string
#include <sys/stat.h>	// chmod()
#include <unistd.h>	// ... unlink()
#include <utility>	// pair<>
#include <vector>

static void grid_start(const char* const prog, const ReadsCorrectionOptions &options, const int i) {
        // create grid script, have grid run it
	ReadsCorrectionOptions new_options(options);
	new_options.job_index = i;
	new_options.grid_options = NULL;
	new_options.grid_options_split = NULL;
	std::string name("m2cns.");
	if (i == -1) {
		name += "split";
		new_options.num_threads = 0;	// flag as grid spin-off
	} else {
       		std::ostringstream x;
		x << i;
		name += x.str();
	}
	std::string script_file(name + ".sh");
	// as we might not have write permission, delete it
	unlink(script_file.c_str());
        std::ofstream out;
        open_fstream(out, script_file.c_str(), std::ios::out);
        //out << "#!/bin/bash\nulimit -c 0\n" << prog << make_options(new_options) << "\n";
        out << "#!/bin/bash\n" << prog << make_options(new_options) << "\n";
	if (!out) {
		std::cerr << "Error writing to " << script_file << "\n";
		exit(1);
	}
        close_fstream(out);
	chmod(script_file.c_str(), S_IRUSR | S_IXUSR);
	std::string cmd(i == -1 && options.grid_options_split ? options.grid_options_split : options.grid_options);
	cmd += " " + name + " " + script_file;
	assert(system(cmd.c_str()) == 0);
}

// pass by value so we can modify list to easily avoid checking previously
// found results files
static void wait_for_files(std::list<std::string> results) {
	const std::list<std::string>::const_iterator end_a(results.end());
	while (!results.empty()) {
		sleep(60);
		std::list<std::string>::iterator a(results.begin());
		while (a != end_a) {
			if (access(a->c_str(), F_OK) == 0) {
				a = results.erase(a);
			} else {
				++a;
			}
		}
	}
}

static void merge_results(const char* const output, const std::list<std::string>& files) {
	std::string out_tmp(output);
	out_tmp += ".tmp";
	std::ofstream out(out_tmp.c_str());
	std::list<std::string>::const_iterator a(files.begin());
	const std::list<std::string>::const_iterator end_a(files.end());
	for (; a != end_a; ++a) {
		std::ifstream in(a->c_str());
		out << in.rdbuf();
		if (!in) {
			std::cerr << "Error reading from " << *a << "\n";
			exit(1);
		} else if (!out) {
			std::cerr << "Error writing to " << out_tmp << "\n";
			exit(1);
		}
	}
	out.close();
	if (rename(out_tmp.c_str(), output) == -1) {
		std::cerr << "Could not rename concatenated output file: " << out_tmp << "\n";
		exit(1);
	}
}

int main(int argc, char** argv) {
	ReadsCorrectionOptions rco;
	if (parse_arguments(argc, argv, rco)) {
		print_usage(argv[0]);
		exit(1);
	} else if (rco.print_usage_info) {
		print_usage(argv[0]);
		exit(0);
	}
	// partition once up front, to allow for grid jobs
	if (rco.job_index == -1) {
		if (access(rco.corrected_reads, F_OK) == 0) {	// full run already done
			return 0;
		} else if (access("partition.done", F_OK) == 0) {	// split already done
		} else if (rco.grid_options || rco.grid_options_split) {
			grid_start(argv[0], rco, -1);
			std::list<std::string> partition_results;
			partition_results.push_back("partition.done");
			wait_for_files(partition_results);
		} else {
			if (rco.input_type == INPUT_TYPE_CAN) {
				// this speeds up candidate starts
				const idx_t n_reads(PackedDB::convert_fasta_to_db(rco.reads, "fasta.db", rco.min_size));
				if (rco.reads_to_correct <= 0 || n_reads < rco.reads_to_correct) {
					rco.reads_to_correct = n_reads;
				}
				// can't reorder reads for memory efficiency - takes
				// too much memory to hold read-read pairings
				partition_candidates(rco.m4, "fasta.db", rco.batch_size, rco.num_partition_files, rco.reads_to_correct);
			} else {
				partition_m4records(rco.m4, rco.min_mapping_ratio - 0.02, rco.batch_size, rco.min_size, rco.num_partition_files);
			}
			if (rco.num_threads == 0) {	// flag saying we're a grid run
				return 0;
			}
		}
	}
	if (rco.grid_options == NULL) {		// single process run
		if (rco.input_type == INPUT_TYPE_CAN) {
			return reads_correction_can(rco);
		} else {
			return reads_correction_m4(rco);
		}
	} else {				// grid run
		// get number of partitions
		std::string idx_file_name;
		generate_partition_index_file_name(rco.m4, idx_file_name);
		std::vector<std::string> partition_file_vec;
		load_partition_files_info(idx_file_name.c_str(), partition_file_vec);
		std::list<std::string> results;
		for (size_t i(0); i != partition_file_vec.size(); ++i) {
			std::ostringstream os;
			os << rco.corrected_reads << "." << i;
			results.push_back(os.str());
			if (access(os.str().c_str(), F_OK) != 0) {
				grid_start(argv[0], rco, i);
				if (rco.grid_start_delay) {
					sleep(rco.grid_start_delay);
				}
			}
		}
		wait_for_files(results);
		merge_results(rco.corrected_reads, results);
	}
	return 0;
}
