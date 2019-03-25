#include "reads_correction_can.h"

#include <cstring>
#include <string>
#include <sstream>

#include "MECAT_AlnGraphBoost.H"
#include "mecat_correction.h"
#include "overlaps_partition.h"
#include "overlaps_store.h"

static void* reads_correction_func_can(void* const arg) {
	ConsensusThreadData& data(*(static_cast<ConsensusThreadData*>(arg)));
	const int tid(data.get_thread_id());
	ConsensusPerThreadData& pdata(data.data[tid]);
	const ExtensionCandidate* const candidates(pdata.candidates);
	const idx_t min_size(ceil(data.rco.min_size * 0.95));
	const int tech_is_pacbio(data.rco.tech == TECH_PACBIO ? 1 : 0);
	idx_t i(pdata.next_candidate);
	while (i < pdata.num_candidates) {
		const idx_t start(i);
		const idx_t sid(candidates[start].sid);
		for (++i; i < pdata.num_candidates && candidates[i].sid == sid; ++i) { }
		if (i - start < data.rco.min_cov || candidates[start].ssize < min_size) {
			continue;
		}
		if (tech_is_pacbio) {
			ns_meap_cns::consensus_one_read_can_pacbio(data, pdata, sid, start, i);
		} else {
			ns_meap_cns::consensus_one_read_can_nanopore(data, pdata, sid, start, i);
		}
		if (pdata.cns_results.size() >= MAX_CNS_RESULTS) {
			data.write_buffer(tid, i);
		}
	}
	data.write_buffer(tid, i);
	return NULL;
}

// load and sort partition data, assign to threads, start threads
static void consensus_one_partition_can(const char* const m4_file_name, const idx_t min_read_id, const idx_t max_read_id, ConsensusThreadData& data) {
	idx_t nec;
	ExtensionCandidate* const ec_list(load_partition_data<ExtensionCandidate>(m4_file_name, nec));
	std::sort(ec_list, ec_list + nec, CmpExtensionCandidateBySidAndScore());
	allocate_ecs(data, ec_list, nec);
	pthread_t thread_ids[data.rco.num_threads];
	for (int i(0); i < data.rco.num_threads; ++i) {
		pthread_create(&thread_ids[i], NULL, reads_correction_func_can, static_cast<void*>(&data));
	}
	for (int i(0); i < data.rco.num_threads; ++i) {
		pthread_join(thread_ids[i], NULL);
	}
	delete[] ec_list;
}

// set up output file and thread data, handle restart if needed
static int reads_correction_can_p(ReadsCorrectionOptions& rco, std::vector<PartitionFileInfo> &partition_file_vec, PackedDB& reads) {
	const PartitionFileInfo& p(partition_file_vec[rco.job_index]);
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
	sprintf(process_info, "processing %s", p.file_name.c_str());
	DynamicTimer dtimer(process_info);
	consensus_one_partition_can(p.file_name.c_str(), p.min_seq_id, p.max_seq_id, data);
	close_fstream(out);
	assert(rename(results_file_tmp.c_str(), results_file.c_str()) == 0);
	return 0;
}

// load reads, partition file index
int reads_correction_can(ReadsCorrectionOptions& rco) {
	std::string idx_file_name;
	generate_partition_index_file_name(rco.m4, idx_file_name);
	std::vector<PartitionFileInfo> partition_file_vec;
	load_partition_files_info(idx_file_name.c_str(), partition_file_vec);
	PackedDB reads;
	reads.load_fasta_db(rco.reads);
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
		for (; rco.job_index < job_end; ++rco.job_index) {
			const PartitionFileInfo& p(partition_file_vec[rco.job_index]);
			sprintf(process_info, "processing %s", p.file_name.c_str());
			DynamicTimer dtimer(process_info);
			consensus_one_partition_can(p.file_name.c_str(), p.min_seq_id, p.max_seq_id, data);
		}
		close_fstream(out);
		assert(rename(tmp_file.c_str(), rco.corrected_reads) == 0);
	}
	return 0;
}
