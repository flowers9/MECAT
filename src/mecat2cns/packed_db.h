#ifndef PACKED_DB_H
#define PACKED_DB_H

#include <fstream>	// ifstream
#include <set>		// set<>
#include <string>	// string
#include <utility>	// pair<>
#include <vector>	// vector<>

#include "../common/defs.h"
#include "../common/sequence.h"		// Sequence
#include "../common/alignment.h"	// ExtensionCandidateCompressed

class PackedDB {
    private:
	struct SeqIndex {
		off_t file_offset;
		idx_t memory_offset, size;
		explicit SeqIndex() { }
		explicit SeqIndex(const off_t i, const idx_t j, const idx_t k) : file_offset(i), memory_offset(j), size(k) { }
	};
    public:
	explicit PackedDB() : pac(0), db_size(0), max_db_size(0), max_read_size_(0) { }
	~PackedDB() {
		delete[] pac;
	}
	// returns number of reads
	static size_t convert_fasta_to_db(const std::string& fasta, const std::string& output_prefix, idx_t min_size);
	static void read_sizes(const std::string& output_prefix, std::vector<idx_t>& sizes);
	// only call one of load_fasta_db and open_db exactly once
	void load_fasta_db(const char* fasta);
	// opens data file, reads in index file
	void open_db(const std::string& filename, idx_t memory_footprint);
	// returns number of candidates that can be processed
	idx_t load_reads(const ExtensionCandidateCompressed* ec_list, idx_t nec);
	void GetSequence(const idx_t id, const bool forward, std::string& seq) const {
		const SeqIndex &si(seq_idx[id]);
		seq.resize(si.size);
		if (forward) {
			const idx_t offset(si.memory_offset);
			for (idx_t i(0); i < si.size; ++i) {
				seq[i] = get_char(offset + i);
			}
		} else {
			const idx_t offset(si.memory_offset + si.size - 1);
			for (idx_t i(0); i < si.size; ++i) {
				seq[i] = 3 - get_char(offset - i);
			}
		}
	}
	idx_t num_reads() const {
		return seq_idx.size();
	}
	idx_t read_size(const idx_t read_id) const {
		return seq_idx[read_id].size;
	}
	idx_t max_read_size() const {
		return max_read_size_;
	}
    private:
	static void set_char(std::vector<uint1>& p, const idx_t idx, const u1_t c) {
		p[idx >> 2] |= c << ((~idx & 3) << 1);
	}
	void set_char(const idx_t idx, const u1_t c) {
		// use ~x instead of 3 - x for speed, since we have to & 3 anyway
		pac[idx >> 2] |= c << ((~idx & 3) << 1);
	}
	u1_t get_char(const idx_t idx) const {
		return pac[idx >> 2] >> ((~idx & 3) << 1) & 3;
	}
	void add_one_seq(const Sequence& seq);
	void destroy() {
		delete[] pac;
		seq_idx.clear();
		max_db_size = db_size = 0;
	}
    private:
	u1_t* pac;
	idx_t db_size, max_db_size, max_read_size_;
	std::vector<SeqIndex> seq_idx;
	std::ifstream pstream;
};

#endif // PACKED_DB_H
