#include "overlaps_partition.h"

#include <fstream>
#include <limits>
#include <set>
#include <vector>

#include "overlaps_store.h"
#include "reads_correction_aux.h"

#define error_and_exit(msg) { std::cerr << msg << "\n"; abort(); }

inline static bool query_is_contained(const M4Record& m4, const double min_cov_ratio) {
	return m4qend(m4) - m4qoff(m4) >= m4qsize(m4) * min_cov_ratio;
}

inline static bool subject_is_contained(const M4Record& m4, const double min_cov_ratio) {
	return m4send(m4) - m4soff(m4) >= m4ssize(m4) * min_cov_ratio;
}

inline static bool check_m4record_mapping_range(const M4Record& m4, const double min_cov_ratio) {
	return query_is_contained(m4, min_cov_ratio) || subject_is_contained(m4, min_cov_ratio);
}

static idx_t get_qualified_m4record_counts(const char* const m4_file_name, const double min_cov_ratio) {
	std::ifstream in;
	open_fstream(in, m4_file_name, std::ios::in);
	idx_t num_reads(-1), num_records(0), num_qualified_records(0);
	M4Record m4;
	m4qext(m4) = m4sext(m4) = INVALID_IDX;
	while (in >> m4) {
		if (m4qext(m4) == INVALID_IDX || m4sext(m4) == INVALID_IDX) {
			ERROR("no gapped start position is provided, please make sure that you have run \'meap_pairwise\' with option \'-g 1\'");
		}
		++num_records;
		if (check_m4record_mapping_range(m4, min_cov_ratio)) {
			++num_qualified_records;
		}
		num_reads = std::max(num_reads, m4qid(m4));
		num_reads = std::max(num_reads, m4sid(m4));
	}
	close_fstream(in);
	LOG(stderr, "there are %ld overlaps, %ld are qualified.", num_records, num_qualified_records);
	return num_reads + 1;
}

// not in use at the moment
#if 0
static void get_repeat_reads(const char* const m4_file_name, const double min_cov_ratio, const idx_t num_reads, std::set<idx_t>& repeat_reads) {
	const int max_contained(100);
	// used to increment to a max of max_contained
	char count_table[max_contained + 1];
	for (int i(0); i < max_contained; ++i) {
		count_table[i] = i + 1;
	}
	count_table[max_contained] = max_contained;
	std::vector<char> counts(num_reads, 0);
	std::ifstream in;
	open_fstream(in, m4_file_name, std::ios::in);
	M4Record m4;
	m4qext(m4) = m4sext(m4) = INVALID_IDX;
	while (in >> m4) {
		if (query_is_contained(m4, min_cov_ratio)) {
			const idx_t qid = m4qid(m4);
			// increments count up to max_contained
			counts[qid] = count_table[static_cast<int>(counts[qid])];
		}
		if (subject_is_contained(m4, min_cov_ratio)) {
			const idx_t sid = m4sid(m4);
			// increments count up to max_contained
			counts[sid] = count_table[static_cast<int>(counts[sid])];
		}
	}
	close_fstream(in);
	for (idx_t i(0); i < num_reads; ++i) {
		if (counts[i] == max_contained) {
			std::cerr << "repeat read " << i << "\n";
			repeat_reads.insert(i);
		}
	}
	LOG(stderr, "number of repeat reads: %lu", repeat_reads.size());
}
#endif

void generate_partition_index_file_name(const char* const m4_file_name, std::string& ret) {
	ret = m4_file_name;
	ret += ".partition_files";
}

void generate_partition_file_name(const char* const m4_file_name, const idx_t part, std::string& ret) {
	ret = m4_file_name;
	ret += ".part";
	std::ostringstream os;
	os << part;
	ret += os.str();
}

static idx_t get_num_reads(const char* const candidates_file) {
	std::ifstream in;
	open_fstream(in, candidates_file, std::ios::in);
	ExtensionCandidate ec;
	int max_id(-1);
	while (in >> ec) {
		max_id = std::max(ec.qid, max_id);
		max_id = std::max(ec.sid, max_id);
	}
	close_fstream(in);
	return max_id + 1;
}

static void normalise_candidate(const ExtensionCandidate& src, ExtensionCandidate& dst, const bool subject_is_target) {
	if (subject_is_target) {
		dst = src;
	} else {
		dst.qdir = src.sdir;
		dst.qid = src.sid;
		dst.qext = src.sext;
		dst.qsize = src.ssize;
		dst.qoff = src.soff;
		dst.qend = src.send;
		dst.sdir = src.qdir;
		dst.sid = src.qid;
		dst.sext = src.qext;
		dst.ssize = src.qsize;
		dst.soff = src.qoff;
		dst.send = src.qend;
		dst.score = src.score;
	}
	if (dst.sdir == REV) {
		dst.qdir = REVERSE_STRAND(dst.qdir);
		dst.sdir = REVERSE_STRAND(dst.sdir);
	}
}

void partition_candidates(const char* input, const idx_t batch_size, const int min_read_size, const int num_files, idx_t num_reads) {
	DynamicTimer dtimer(__func__);
	PartitionResultsWriter<ExtensionCandidate> prw(num_files);
	idx_t i(0);
	off_t input_pos;
	int is_restart(prw.restart(input, generate_partition_file_name, "partition.done", i, input_pos));
	if (!is_restart) {
		prw.num_reads = num_reads ? num_reads : get_num_reads(input);
	}
	const idx_t num_batches((prw.num_reads + batch_size - 1) / batch_size);
	std::string idx_file_name;
	generate_partition_index_file_name(input, idx_file_name);
	std::ofstream idx_file;
	open_fstream(idx_file, idx_file_name.c_str(), std::ios::out);
	ExtensionCandidate ec, nec;
	// and here we go through the input file num_batches times,
	// being limited by the number of open output files we can have
	for (; i < num_batches; i += prw.kNumFiles) {
		const idx_t sfid(i);
		const idx_t efid(std::min(sfid + prw.kNumFiles, num_batches));
		const int nf(efid - sfid);
		const idx_t L(batch_size * sfid);
		const idx_t R(efid < num_batches ? batch_size * efid : prw.num_reads);
		std::ifstream in;
		open_fstream(in, input, std::ios::in);
		if (is_restart) {
			if (!in.seekg(input_pos)) {
				ERROR("Input seek failed while restoring checkpoint: %s", input);
			}
			is_restart = 0;
		} else {
			prw.OpenFiles(sfid, efid, input, generate_partition_file_name, "partition.done");
		}
		while (in >> ec) {
			if (ec.qsize < min_read_size || ec.ssize < min_read_size) {
				continue;
			}
			// not set by >>
			ec.qoff = ec.soff = ec.qend = ec.send = 0;
			if (L <= ec.qid && ec.qid < R) {
				normalise_candidate(ec, nec, false);
				if (prw.WriteOneResult((ec.qid - L) / batch_size, ec.qid, nec)) {
					prw.checkpoint(in.tellg());
				}
			}
			if (L <= ec.sid && ec.sid < R) {
				normalise_candidate(ec, nec, true);
				if (prw.WriteOneResult((ec.sid - L) / batch_size, ec.sid, nec)) {
					prw.checkpoint(in.tellg());
				}
			}
		}
		for (int k(0); k < nf; ++k) {
			fprintf(stderr, "%s contains %ld overlaps\n", prw.file_names[k].c_str(), prw.counts[k]);
			if (prw.counts[k] != 0) {
				idx_file << prw.file_names[k] << "\n";
			}
		}
		prw.CloseFiles();
	}
	close_fstream(idx_file);
	prw.finalize();
}

void partition_m4records(const char* const m4_file_name, const double min_cov_ratio, const idx_t batch_size, const int min_read_size, const int num_files) {
	DynamicTimer dtimer(__func__);
	idx_t num_reads(get_qualified_m4record_counts(m4_file_name, min_cov_ratio));
	std::set<idx_t> repeat_reads;
	//get_repeat_reads(m4_file_name, min_cov_ratio, num_reads, repeat_reads);
	const idx_t num_batches((num_reads + batch_size - 1) / batch_size);
	std::string idx_file_name;
	generate_partition_index_file_name(m4_file_name, idx_file_name);
	std::ofstream idx_file;
	open_fstream(idx_file, idx_file_name.c_str(), std::ios::out);
	M4Record m4, nm4;
	ExtensionCandidate ec;
	PartitionResultsWriter<ExtensionCandidate> prw(num_files);
	for (idx_t i(0); i < num_batches; i += prw.kNumFiles) {
		const idx_t sfid(i);
		const idx_t efid(std::min(sfid + prw.kNumFiles, num_batches));
		const int nf(efid - sfid);
		const idx_t L(batch_size * sfid);
		const idx_t R(efid < num_batches ? batch_size * efid : prw.num_reads);
		std::ifstream in;
		open_fstream(in, m4_file_name, std::ios::in);
		prw.OpenFiles(sfid, efid, m4_file_name, generate_partition_file_name, "partition.done");
		while (in >> m4) {
			if (m4qsize(m4) < min_read_size || m4ssize(m4) < min_read_size) {
				continue;
			} else if (!check_m4record_mapping_range(m4, min_cov_ratio)) {
				continue;
			} else if (repeat_reads.find(m4qid(m4)) != repeat_reads.end() || repeat_reads.find(m4sid(m4)) != repeat_reads.end()) {
				continue;
			}
			if (m4qid(m4) >= L && m4qid(m4) < R) {
				normalize_m4record(m4, false, nm4);
				m4_to_candidate(nm4, ec);
				prw.WriteOneResult((m4qid(m4) - L) / batch_size, m4qid(m4), ec);
			}
			if (m4sid(m4) >= L && m4sid(m4) < R) {
				normalize_m4record(m4, true, nm4);
				m4_to_candidate(nm4, ec);
				prw.WriteOneResult((m4sid(m4) - L) / batch_size, m4sid(m4), ec);
			}
		}
		for (int k(0); k < nf; ++k) {
			fprintf(stderr, "%s contains %ld overlaps\n", prw.file_names[k].c_str(), prw.counts[k]);
			if (prw.counts[k] != 0) {
				idx_file << prw.file_names[k] << "\n";
			}
		}
		prw.CloseFiles();
	}
	close_fstream(idx_file);
}

void load_partition_files_info(const char* const idx_file_name, std::vector<std::string>& file_info_vec) {
	file_info_vec.clear();
	std::ifstream in;
	open_fstream(in, idx_file_name, std::ios::in);
	std::string line;
	while (in >> line) {
		file_info_vec.push_back(line);
	}
	close_fstream(in);
}
