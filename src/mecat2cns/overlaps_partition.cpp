#include "overlaps_partition.h"

#include <fstream>
#include <limits>
#include <set>
#include <sys/stat.h>	// stat(), struct stat
#include <vector>
#include <utility>	// make_pair(), pair<>

#include "overlaps_store.h"
#include "reads_correction_aux.h"
#include "packed_db.h"	// PackedDB

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

void partition_candidates(const std::string& input, const std::string& pac_prefix, const size_t batch_size, const int num_files, const idx_t num_reads) {
	DynamicTimer dtimer(__func__);
	struct stat buf;
	if (stat(input.c_str(), &buf) == -1) {
		ERROR("Could not get file size: %s", input.c_str());
	}
	// each candidate line takes on average 44 characters, but go with 32;
	// each one produces two candidates (forward and reverse)
	const idx_t num_batches(buf.st_size / 32 * 2 * sizeof(ExtensionCandidateCompressed) / batch_size);
	// separate them by read id on this pass, as we don't know how
	// many candidates each read is part of; we're assuming an even
	// distribution on average
	const idx_t reads_per_batch((num_reads + num_batches - 1) / num_batches);
	std::vector<idx_t> read_sizes;
	PackedDB::read_sizes(pac_prefix, read_sizes);
	PartitionResultsWriter<ExtensionCandidateCompressed> prw(num_files);
	idx_t i(0);
	off_t input_pos;
	int is_restart(prw.restart(input, generate_partition_file_name, "partition.done", i, input_pos));
	std::string idx_file_name;
	generate_partition_index_file_name(input.c_str(), idx_file_name);
	std::ofstream idx_file;
	open_fstream(idx_file, idx_file_name.c_str(), std::ios::out);
	ExtensionCandidate ec;
	ExtensionCandidateCompressed nec;
	// and here we go through the input file num_batches times,
	// being limited by the number of open output files we can have
	for (; i < num_batches; i += prw.kNumFiles) {
		const idx_t sfid(i);
		const idx_t efid(std::min(sfid + prw.kNumFiles, num_batches));
		const int nf(efid - sfid);
		const idx_t L(sfid * reads_per_batch);
		const idx_t R(efid < num_batches ? efid * reads_per_batch : num_reads);
		std::ifstream in;
		open_fstream(in, input.c_str(), std::ios::in);
		if (is_restart) {
			if (!in.seekg(input_pos)) {
				ERROR("Input seek failed while restoring checkpoint: %s", input.c_str());
			}
			is_restart = 0;
		} else {
			prw.OpenFiles(sfid, efid, input, generate_partition_file_name, "partition.done");
		}
		while (in >> ec) {
			if (!read_sizes[ec.sid] || !read_sizes[ec.qid]) {
				continue;
			}
			// make sure these match before tossing them
			r_assert(ec.ssize == read_sizes[ec.sid] && ec.qsize == read_sizes[ec.qid]);
			// as variable sizes may be different, make sure the ones
			// we read in can be safely stored
			r_assert(ec.sid <= ExtensionCandidateCompressed::max_value);
			r_assert(ec.qid <= ExtensionCandidateCompressed::max_value);
			r_assert(ec.sext <= ExtensionCandidateCompressed::max_value);
			// qext is one bit smaller than the others
			r_assert(ec.qext <= ExtensionCandidateCompressed::max_qext);
			r_assert(ec.score <= ExtensionCandidateCompressed::max_value);
			if (L <= ec.sid && ec.sid < R) {
				nec.set(ec);
				if (prw.WriteOneResult((nec.sid - L) / reads_per_batch, nec.sid, nec)) {
					prw.checkpoint(in.tellg());
				}
			}
			if (L <= ec.qid && ec.qid < R) {
				nec.set_swap(ec);
				if (prw.WriteOneResult((nec.sid - L) / reads_per_batch, nec.sid, nec)) {
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

// assign reads to files in read order so that each file
// has about the same number of candidates

static void allocate_reads_to_files(const idx_t num_batches, const std::vector<int>& align_counts, std::vector<int>& read_to_file) {
	idx_t total_aligns(0);
	std::vector<int>::const_iterator a(align_counts.begin());
	const std::vector<int>::const_iterator end_a(align_counts.end());
	for (; a != end_a; ++a) {
		total_aligns += *a;
	}
	idx_t read_id(0);
	for (int batch(0); batch != num_batches; ++batch) {
		idx_t align_count(0);
		// drop fractions here, as we'll overcount some below
		const idx_t want_aligns(total_aligns / (num_batches - batch));
		for (; align_count < want_aligns; ++read_id) {
			// skip over reads we're not correcting, as they might
			// actually have alignments (which we want to ignore)
			if (align_counts[read_id]) {
				align_count += align_counts[read_id];
				read_to_file[read_id] = batch;
			}
		}
		total_aligns -= align_count;
	}
}

void partition_candidates_reorder(const std::string& input, const size_t batch_size, const int num_files, const idx_t num_reads, const std::vector<idx_t>& read_order, const std::vector<int>& align_counts) {
	DynamicTimer dtimer(__func__);
	struct stat buf;
	if (stat(input.c_str(), &buf) == -1) {
		ERROR("Could not get file size: %s", input.c_str());
	}
	// each candidate line takes approximately 44 characters, but go with 32;
	// each one produces two candidates (forward and reverse)
	const int num_batches(buf.st_size / 32 * 2 * sizeof(ExtensionCandidateCompressed) / batch_size);
	PartitionResultsWriter<ExtensionCandidateCompressed> prw(num_files);
	idx_t i(0);
	off_t input_pos;
	int is_restart(prw.restart(input, generate_partition_file_name, "partition.done", i, input_pos));
	std::vector<int> read_to_file(num_reads, -1);
	allocate_reads_to_files(num_batches, align_counts, read_to_file);
	std::string idx_file_name;
	generate_partition_index_file_name(input.c_str(), idx_file_name);
	std::ofstream idx_file;
	open_fstream(idx_file, idx_file_name.c_str(), std::ios::out);
	ExtensionCandidate ec;
	ExtensionCandidateCompressed nec;
	// and here we go through the input file num_batches times,
	// being limited by the number of open output files we can have
	for (; i < num_batches; i += prw.kNumFiles) {
		const int sfid(i);
		const int efid(std::min(sfid + prw.kNumFiles, num_batches));
		std::ifstream in;
		open_fstream(in, input.c_str(), std::ios::in);
		if (is_restart) {
			if (!in.seekg(input_pos)) {
				ERROR("Input seek failed while restoring checkpoint: %s", input.c_str());
			}
			is_restart = 0;
		} else {
			prw.OpenFiles(sfid, efid, input, generate_partition_file_name, "partition.done");
		}
		while (in >> ec) {
			ec.qid = read_order[ec.qid];	// convert to new read order
			ec.sid = read_order[ec.sid];
			if (ec.qid == -1 || ec.sid == -1) {
				continue;
			}
			r_assert(ec.sid <= ExtensionCandidateCompressed::max_value);
			r_assert(ec.qid <= ExtensionCandidateCompressed::max_value);
			r_assert(ec.sext <= ExtensionCandidateCompressed::max_value);
			r_assert(ec.qext <= ExtensionCandidateCompressed::max_qext);
			r_assert(ec.score <= ExtensionCandidateCompressed::max_value);
			const int sfile(read_to_file[ec.sid]);
			if (sfid <= sfile && sfile < efid) {
				nec.set(ec);
				if (prw.WriteOneResult(sfile, nec.sid, nec)) {
					prw.checkpoint(in.tellg());
				}
			}
			const int qfile(read_to_file[ec.qid]);
			if (sfid <= qfile && qfile < efid) {
				nec.set_swap(ec);
				if (prw.WriteOneResult(qfile, nec.sid, nec)) {
					prw.checkpoint(in.tellg());
				}
			}
		}
		const int nf(efid - sfid);
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

void partition_m4records(const char* const m4_file_name, const double min_cov_ratio, const size_t batch_size, const int min_read_size, const int num_files) {
	DynamicTimer dtimer(__func__);
	idx_t num_reads(get_qualified_m4record_counts(m4_file_name, min_cov_ratio));
	std::set<idx_t> repeat_reads;
	//get_repeat_reads(m4_file_name, min_cov_ratio, num_reads, repeat_reads);
	struct stat buf;
	if (stat(m4_file_name, &buf) == -1) {
		ERROR("Could not get file size: %s", m4_file_name);
	}
	const idx_t num_batches(buf.st_size / 32 * 2 * sizeof(ExtensionCandidate) / batch_size);
	const idx_t reads_per_batch((num_reads + num_batches - 1) / num_batches);
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
		const idx_t L(reads_per_batch * sfid);
		const idx_t R(efid < num_batches ? reads_per_batch * efid : num_reads);
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
				prw.WriteOneResult((m4qid(m4) - L) / reads_per_batch , m4qid(m4), ec);
			}
			if (m4sid(m4) >= L && m4sid(m4) < R) {
				normalize_m4record(m4, true, nm4);
				m4_to_candidate(nm4, ec);
				prw.WriteOneResult((m4sid(m4) - L) / reads_per_batch , m4sid(m4), ec);
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

// XXX - convert this to using a file split using PartitionResultsWriter, as it's
// just way too big to hold in memory for 1TB of reads, plus it takes forever
// to go through the 10TB alignment file, so that just has to be checkpointed,
// plus the whole point is that we're memory limited when we go this way, so why
// gamble on having enough memory to hold the array?

static void generate_new_read_order(const std::string& input, const idx_t num_reads, const int min_cov, const std::vector<idx_t>& old_read_sizes, std::vector<idx_t>& new_order, std::vector<int>& align_counts) {
	// first read in candidates and find all read-read pairings
	std::vector<std::pair<idx_t, idx_t> > aligns;
	ExtensionCandidate ec;
	std::ifstream in;
	open_fstream(in, input.c_str(), std::ios::in);
	while (in >> ec) {
		// screen out small reads, aligns between reads we don't care about;
		// use old_read_sizes as we've already screened out the small reads
		if (old_read_sizes[ec.sid] && old_read_sizes[ec.qid] && (ec.qid < num_reads || ec.sid < num_reads)) {
			r_assert(old_read_sizes[ec.sid] == ec.ssize && old_read_sizes[ec.qid] == ec.qsize);
			if (ec.qid < num_reads) {
				aligns.push_back(std::make_pair(ec.qid, ec.sid));
			}
			if (ec.sid < num_reads) {
				aligns.push_back(std::make_pair(ec.sid, ec.qid));
			}
		}
	}
	close_fstream(in);
	std::sort(aligns.begin(), aligns.end());
	// generate index into aligns
	align_counts.assign(num_reads, 0);
	std::vector<idx_t> aligns_index(num_reads, -1);
	const idx_t end_i(aligns.size());
	for (idx_t i(0); i != end_i;) {
		const idx_t start(i);
		const idx_t read_id(aligns[i].first);
		for (++i; i != end_i && aligns[i].first == read_id; ++i) { }
		if (i - start >= min_cov) {
			aligns_index[read_id] = start;
			align_counts[read_id] = i - start;
		}
	}
	// generate the new read order;
	// make sure to make space for all reads that
	// could be used, not just those we're correcting
	std::vector<char> used(old_read_sizes.size(), 0);
	new_order.reserve(old_read_sizes.size());
	size_t next_search(0);
	for (idx_t next_unused(0);;) {
		// skip over used reads, reads with too few alignments
		for (; next_unused != num_reads && (used[next_unused] || aligns_index[next_unused] == -1); ++next_unused) { }
		if (next_unused == num_reads) {
			break;
		}
		used[next_unused] = 1;
		new_order.push_back(next_unused);
		// add all reads aligned to, and aligned to those, and so on
		for (; next_search != new_order.size(); ++next_search) {
			const idx_t sid(new_order[next_search]);
			// we only use aligns to first num_reads reads
			if (sid < num_reads) {
				idx_t i(aligns_index[sid]);
				if (i != -1) {
					for (; i != end_i && aligns[i].first == sid; ++i) {
						const idx_t qid(aligns[i].second);
						if (!used[qid]) {
							used[qid] = 1;
							new_order.push_back(qid);
						}
					}
				}
			}
		}
	}
}

// create a new order for reads to speed up pulling in reads for alignment
// processing (for conditions where memory is not sufficient to hold all
// reads in memory); generates a file with the old -> new read ordering,
// also generates an index for a fasta db file with the reads in the new
// order (to allow easy conversion of the original fasta to a fasta db
// with the new read ordering)

void make_read_sort_order(const std::string& input, const std::string& old_pac_prefix, const std::string& sort_file_name, const std::string& pac_prefix, const idx_t num_reads, const int min_cov, std::vector<idx_t>& read_order, std::vector<std::pair<idx_t, idx_t> >& read_index, std::vector<int>& align_counts) {
	// if file already exists, just read it in
	if (access(sort_file_name.c_str(), F_OK) == 0) {
		struct stat buf;
		if (stat(sort_file_name.c_str(), &buf) == -1) {
			ERROR("Could not stat read reorder file: %s", sort_file_name.c_str());
		}
		read_order.resize(buf.st_size / sizeof(idx_t));
		std::ifstream in;
		open_fstream(in, sort_file_name.c_str(), std::ios::in);
		if (!in.read((char *)(&read_order[0]), buf.st_size)) {
			ERROR("Error reading read reorder file: %s", sort_file_name.c_str());
		}
		close_fstream(in);
		PackedDB::read_index(pac_prefix, read_index);
		return;
	}
	DynamicTimer dtimer(__func__);
	std::vector<idx_t> old_read_sizes;	// pre-reorder sizes
	PackedDB::read_sizes(old_pac_prefix, old_read_sizes);
	std::vector<idx_t> new_order;		// [new_read_id] = old_read_id
	std::vector<int> presort_align_counts;
	generate_new_read_order(input, num_reads, min_cov, old_read_sizes, new_order, presort_align_counts);
	// now reverse new_order into read_order;
	// also, generate index for reordered read database
	read_order.assign(old_read_sizes.size(), -1);
	read_index.resize(new_order.size());
	align_counts.assign(new_order.size(), 0);
	idx_t total_size(0);
	const idx_t end_new_rid(new_order.size());
	for (idx_t new_rid(0); new_rid != end_new_rid; ++new_rid) {
		const idx_t old_rid(new_order[new_rid]);
		read_order[old_rid] = new_rid;
		std::pair<idx_t, idx_t>& b(read_index[new_rid]);
		b.first = total_size;
		b.second = old_read_sizes[old_rid];
		total_size += (b.second + 3) / 4;
		if (old_rid < num_reads) {
			align_counts[new_rid] = presort_align_counts[old_rid];
		}
	}
	PackedDB::create_index(pac_prefix, read_index);
	const std::string sort_file_name_tmp(sort_file_name + ".tmp");
	std::ofstream out;
	open_fstream(out, sort_file_name_tmp.c_str(), std::ios::out | std::ios::binary);
	if (!out.write((char *)(&read_order[0]), sizeof(idx_t) * read_order.size())) {
		ERROR("Error writing to read reorder file: %s", sort_file_name_tmp.c_str());
	}
	close_fstream(out);
	if (rename(sort_file_name_tmp.c_str(), sort_file_name.c_str()) == -1) {
		ERROR("Could not rename read reorder file: %s", sort_file_name_tmp.c_str());
	}
}
