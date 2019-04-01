#ifndef PACKED_DB_H
#define PACKED_DB_H

#include <fstream>	// ifstream
#include <set>		// set<>
#include <string>	// string

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
	static std::string generate_pac_name(const std::string& prefix) {
		return prefix + ".pac";
	}
	static std::string generate_idx_name(const std::string& prefix) {
		return prefix + ".idx";
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

#if 0
	// these routines aren't used anywhere
	void clear() {
		seq_idx.clear();
		db_size = 0;
		memset(pac, 0, (max_db_size + 3) / 4);
	}
	void reserve(const idx_t size) {
		destroy();
		max_db_size = size;
		safe_calloc(pac, u1_t, (max_db_size + 3) / 4);
	}
	void get_sequence(const idx_t from, const idx_t to, const bool forward, char* const seq) const {
		if (forward) {
			idx_t idx(0);
			for (idx_t i(from); i < to; ++i, ++idx) {
				seq[idx] = get_char(pac, i);
			}
		} else {
			idx_t idx(to - from - 1);
			for (idx_t i(from); i < to; ++i, --idx) {
				seq[idx] = 3 - get_char(pac, i);
			}
		}
	}
	void get_sequence(const idx_t rid, const bool forward, char* const seq) const {
		const SeqIndex& a(seq_idx[rid]);
		get_sequence(a.memory_offset, a.memory_offset + a.size, forward, seq);
	}
	void get_sequence(const idx_t rid, const idx_t from, const idx_t to, const bool forward, char* const seq) const {
		const idx_t offset(seq_idx[rid].memory_offset);
		get_sequence(offset + from, offset + to, forward, seq);
	}
	static void decode_sequence(char* const seq, const idx_t seq_size) {
		for (idx_t i(0); i < seq_size; ++i) {
			const u1_t c(seq[i]);
			r_assert(c < 4);	// c is unsigned, so always >= 0
			seq[i] = "ACGT"[c];
		}
	}
	void get_decode_sequence(const idx_t rid, const idx_t from, const idx_t to, const bool forward, char* const seq) const {
		get_sequence(rid, from, to, forward, seq);
		decode_sequence(seq, to - from);
	}
	void set_char(const idx_t idx, const u1_t c) {
		set_char(pac, idx, c);
	}
	u1_t get_char(const idx_t idx) const {
		return get_char(pac, idx);
	}
	idx_t size() const {
		return db_size;
	}
	idx_t num_seqs() const {
		return seq_idx.size();
	}
	idx_t seq_offset(const idx_t rid) const {
		return seq_idx[rid].memory_offset;
	}
	idx_t seq_size(const idx_t rid) const {
		return seq_idx[rid].size;
	}
	void add_one_seq(const char* seq, const idx_t size);
	idx_t offset_to_rid(const idx_t offset) const;
	static void dump_pac(const u1_t* p, idx_t size, const char* path);
	static u1_t* load_pac(const char* path, idx_t& size);
	static void dump_idx(const PODArray<SeqIndex>& idx_list, const char* path);
	static void load_idx(const char* path, PODArray<SeqIndex>& idx_list);
	void dump_packed_db(const char* path) const;
	void load_packed_db(const char* path);
	static void pack_fasta_db(const char* fasta, const char* output_prefix, const idx_t min_size);
#endif

};

#endif // PACKED_DB_H
