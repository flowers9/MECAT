#include "reads_correction_can.h"

#include <cstring>
#include <string>	// string
#include <sstream>

#include "MECAT_AlnGraphBoost.H"
#include "mecat_correction.h"
#include "overlaps_partition.h"
#include "overlaps_store.h"

static void* reads_correction_func_can(void* const arg) {
	ConsensusThreadData& data(*(static_cast<ConsensusThreadData*>(arg)));
	const int tid(data.get_thread_id());
	ConsensusPerThreadData& pdata(data.data[tid]);
	const ExtensionCandidateCompressed* const candidates((ExtensionCandidateCompressed*)pdata.candidates);
	idx_t i(pdata.next_candidate);
	if (data.rco.tech == TECH_PACBIO) {
		while (i != pdata.num_candidates) {
			const idx_t start(i);
			const idx_t sid(candidates[start].sid);
			for (++i; i != pdata.num_candidates && candidates[i].sid == sid; ++i) { }
			// still have to check this here, as it's not always checked earlier
			if (i - start < data.rco.min_cov) {
				continue;
			}
			ns_meap_cns::consensus_one_read_can_pacbio(data, pdata, sid, start, i);
			if (pdata.cns_results.size() >= MAX_CNS_RESULTS) {
				data.write_buffer(tid, i);
			}
		}
	} else {
		while (i != pdata.num_candidates) {
			const idx_t start(i);
			const idx_t sid(candidates[start].sid);
			for (++i; i != pdata.num_candidates && candidates[i].sid == sid; ++i) { }
			if (i - start < data.rco.min_cov) {
				continue;
			}
			ns_meap_cns::consensus_one_read_can_nanopore(data, pdata, sid, start, i);
			if (pdata.cns_results.size() >= MAX_CNS_RESULTS) {
				data.write_buffer(tid, i);
			}
		}
	}
	data.write_buffer(tid, i);
	return NULL;
}

class EC_Index {	// offset into ec_list (and number of ecs) for each read id
    public:
	idx_t offset, count;
	explicit EC_Index() : offset(0), count(0) { }
	explicit EC_Index(idx_t i, idx_t j) : offset(i), count(j) { }
	~EC_Index() { }
};

// reorder list so reads are more concentrated and we can process more
// ecs per pass, with limited read space; also filters list to exclude
// candidates of reads with low coverage

static ExtensionCandidateCompressed* reorder_candidates(ExtensionCandidateCompressed* const ec_list, idx_t& nec, const idx_t num_reads, const int min_cov) {
	idx_t total_ec(0);
	std::vector<EC_Index> index(num_reads);
	// index existing list by sid
	for (idx_t i(0); i != nec;) {
		const idx_t start(i);
		const idx_t sid(ec_list[i].sid);
		for (++i; i != nec && ec_list[i].sid == sid; ++i) { }
		// make sure we have enough coverage
		// (don't need to check size, that happened during partition)
		const idx_t count(i - start);
		if (count >= min_cov) {
			index[sid] = EC_Index(start, count);
			total_ec += count;
		}
	}
	ExtensionCandidateCompressed* new_list(new ExtensionCandidateCompressed[total_ec]);
	// generate the new read order
	std::vector<char> used(num_reads, 0);
	std::vector<idx_t> new_order;
	new_order.reserve(num_reads);
	idx_t next_unused(0);
	size_t next_search(0);
	for (;;) {
		// skip over used reads, reads with no alignments
		for (; next_unused != num_reads && (used[next_unused] || index[next_unused].count == 0); ++next_unused) { }
		if (next_unused == num_reads) {
			break;
		}
		used[next_unused] = 1;
		new_order.push_back(next_unused);
		// add all reads aligned to, and aligned to those, and so on
		for (; next_search != new_order.size(); ++next_search) {
			const idx_t sid(new_order[next_search]);
			const EC_Index& a(index[sid]);
			// only include reads that align to this one if we plan to
			// actually run those alignments (that is, if it has sufficient
			// coverage)
			if (a.count >= min_cov) {
				idx_t i(a.offset);
				const idx_t end_i(i + a.count);
				for (; i != end_i; ++i) {
					const idx_t qid(ec_list[i].qid);
					if (!used[qid]) {
						used[qid] = 1;
						new_order.push_back(qid);
					}
				}
			}
		}
	}
	// copy over the ecs to the new list in the new order
	idx_t pos(0);
	std::vector<idx_t>::const_iterator a(new_order.begin());
	const std::vector<idx_t>::const_iterator end_a(new_order.end());
	for (; a != end_a; ++a) {
		const EC_Index& b(index[*a]);
		if (b.count) {
			memcpy(new_list + pos, ec_list + b.offset, sizeof(ExtensionCandidateCompressed) * b.count);
			pos += b.count;
		}
	}
	delete[] ec_list;
	// only do this at the end, once we know we'll successfully complete
	nec = total_ec;
	return new_list;
}

// load and sort partition data, assign to threads, start threads
static void consensus_one_partition_can(const char* const m4_file_name, ConsensusThreadData& data) {
	idx_t nec;
	ExtensionCandidateCompressed* ec_list(load_partition_data<ExtensionCandidateCompressed>(m4_file_name, nec));
	std::sort(ec_list, ec_list + nec, CmpExtensionCandidateCompressedBySidAndScore());
	// if we're memory limited and we didn't already reorder, do it here;
	// spend some cpu time to reduce number of passes
	if (data.rco.read_buffer_size) {
		// don't die if we run out of memory, just do it the slow way
		try {
			ec_list = reorder_candidates(ec_list, nec, data.reads.num_reads(), data.rco.min_cov);
		} catch (std::bad_alloc &except) { }
	}
	pthread_t thread_ids[data.rco.num_threads];
	while (data.ec_offset != nec) {
		// see how many candidates we can run, given
		// how much read sequence we can load into memory
		// (unless we're not limited, in which case load 'em all)
		const idx_t ecs(data.rco.read_buffer_size ? data.reads.load_reads(ec_list + data.ec_offset, nec - data.ec_offset) : nec);
		allocate_ecs(data, ec_list + data.ec_offset, ecs);
		for (int i(0); i != data.rco.num_threads; ++i) {
			pthread_create(&thread_ids[i], NULL, reads_correction_func_can, static_cast<void*>(&data));
		}
		for (int i(0); i != data.rco.num_threads; ++i) {
			pthread_join(thread_ids[i], NULL);
		}
		// update offset after threading to avoid checkpointing with new value
		data.ec_offset += ecs;
		data.reset_threads();
	}
	delete[] ec_list;
}

// set up output file and thread data, handle restart if needed
static int reads_correction_can_p(ReadsCorrectionOptions& rco, std::vector<std::string> &partition_file_vec, PackedDB& reads) {
	const std::string& p(partition_file_vec[rco.job_index]);
	std::ostringstream os;
	os << rco.corrected_reads << "." << rco.job_index;
	const std::string results_file(os.str());
	const std::string results_file_tmp(results_file + ".tmp");
	std::ofstream out;
	ConsensusThreadData data(rco, reads, out, results_file.c_str());
	off_t output_pos;
	if (data.restart(output_pos)) {					// is a restart
		open_fstream(out, results_file_tmp.c_str(), std::ios::out | std::ios::in);
		if (!out.seekp(output_pos)) {
			ERROR("Restart failed: output file seek failed: %s", results_file_tmp.c_str());
		}
	} else {
		open_fstream(out, results_file_tmp.c_str(), std::ios::out);
	}
	char process_info[1024];
	sprintf(process_info, "processing %s", p.c_str());
	DynamicTimer dtimer(process_info);
	consensus_one_partition_can(p.c_str(), data);
	close_fstream(out);
	assert(rename(results_file_tmp.c_str(), results_file.c_str()) == 0);
	return 0;
}

// load reads, partition file index
int reads_correction_can(ReadsCorrectionOptions& rco) {
	std::string idx_file_name;
	generate_partition_index_file_name(rco.m4, idx_file_name);
	std::vector<std::string> partition_file_vec;
	load_partition_files_info(idx_file_name.c_str(), partition_file_vec);
	PackedDB reads;
	reads.open_db("fasta.db", rco.read_buffer_size);
	if (rco.job_index != -1) {
		return reads_correction_can_p(rco, partition_file_vec, reads);
	} else {
		rco.job_index = 0;	// may get changed by data.restart() below
		std::string tmp_file(rco.corrected_reads);
		tmp_file += ".tmp";
		std::ofstream out;
		ConsensusThreadData data(rco, reads, out, rco.corrected_reads);
		off_t output_pos;
		if (data.restart(output_pos)) {				// is a restart
			open_fstream(out, rco.corrected_reads, std::ios::out | std::ios::in);
			out.seekp(output_pos);
			if (!out.seekp(output_pos)) {
				ERROR("Restart failed: output file seek failed: %s", rco.corrected_reads);
			}
		} else {
			open_fstream(out, rco.corrected_reads, std::ios::out);
		}
		char process_info[1024];
		const int job_end(partition_file_vec.size());
		for (; rco.job_index != job_end; ++rco.job_index) {
			const std::string& p(partition_file_vec[rco.job_index]);
			sprintf(process_info, "processing %s", p.c_str());
			DynamicTimer dtimer(process_info);
			consensus_one_partition_can(p.c_str(), data);
		}
		close_fstream(out);
		assert(rename(tmp_file.c_str(), rco.corrected_reads) == 0);
	}
	return 0;
}
