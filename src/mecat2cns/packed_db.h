#ifndef PACKED_DB_H
#define PACKED_DB_H

#include <fstream>	// ifstream
#include <set>		// set<>
#include <string>	// string
#include <utility>	// pair<>
#include <vector>	// vector<>

#include "../common/defs.h"
#include "../common/sequence.h"		// Sequence
#include "../common/alignment.h"	// ExtensionCandidate

class PackedDB {
    private:
	struct SeqIndex {
		off_t file_offset;
		idx_t memory_offset, size;
	};
    public:
	PackedDB() : pac(NULL), db_size(0), max_db_size(0) { }
	~PackedDB() {
		if (pac) {
			safe_free(pac);
		}
	}
	// only call one of load_fasta_db and open_db exactly once
	void load_fasta_db(const char* fasta);
	static void convert_fasta_to_db(const std::string& fasta, const std::string& output_prefix, idx_t min_size);
	// create a prospective index file for pac file to be written in random order
	static void create_index(const std::string& output_prefix, const std::vector<std::pair<idx_t, idx_t> >& index);
	static void read_index(const std::string& output_prefix, std::vector<std::pair<idx_t, idx_t> >& index);
	static void convert_fasta_to_ordered_db(const std::string& fasta, const std::string& output_prefix, const std::vector<std::pair<idx_t, idx_t> >& index, const std::vector<idx_t>& read_order);
	// opens data file, reads in index file
	void open_db(const std::string& filename, idx_t memory_footprint);
	// returns number of candidates that can be processed
	idx_t load_reads(const ExtensionCandidate* ec_list, idx_t nec);
	void GetSequence(const idx_t id, const bool forward, char* const seq, const idx_t size) const {
		r_assert(size == seq_idx[id].size);
		if (forward) {
			const idx_t offset(seq_idx[id].memory_offset);
			for (idx_t i(0); i < size; ++i) {
				seq[i] = get_char(offset + i);
			}
		} else {
			const idx_t offset(seq_idx[id].memory_offset + size - 1);
			for (idx_t i(0); i < size; ++i) {
				seq[i] = 3 - get_char(offset - i);
			}
		}
	}
	idx_t num_reads() const {
		return seq_idx.size();
	}
    private:
	static void set_char(u1_t* const p, const idx_t idx, const u1_t c) {
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
		if (pac) {
			safe_free(pac);
			pac = NULL;
		}
		seq_idx.clear();
		max_db_size = db_size = 0;
	}
    private:
	u1_t* pac;
	idx_t db_size;
	idx_t max_db_size;
	PODArray<SeqIndex> seq_idx;
	std::ifstream pstream;
};

#endif // PACKED_DB_H
