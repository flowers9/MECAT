#ifndef PACKED_DB_H
#define PACKED_DB_H

#include "defs.h"
#include "sequence.h"

class PackedDB {
    public:
	struct SeqIndex {
		idx_t id, offset, size;
	};
    public:
	PackedDB() : pac(NULL), db_size(0), max_db_size(0) {}
	~PackedDB() {
		destroy();
	}
	void reserve(const idx_t size) {
		destroy();
		max_db_size = size;
		idx_t bytes = max_db_size / 4;
		safe_calloc(pac, u1_t, bytes);
	}
	void GetSequence(const idx_t id, const bool fwd, char* const seq, const idx_t size) const {
		r_assert(size == seq_idx[id].size);
		if (fwd) {
			const idx_t offset(seq_idx[id].offset);
			for (idx_t i(0); i < size; ++i) {
				seq[i] = get_char(offset + i);
			}
		} else {
			const idx_t offset(seq_idx[id].offset + size - 1);
			for (idx_t i(0); i < size; ++i) {
				seq[i] = 3 - get_char(offset - i);
			}
		}
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
		get_sequence(a.offset, a.offset + a.size, forward, seq);
	}
	void get_sequence(const idx_t rid, const idx_t from, const idx_t to, const bool forward, char* const seq) const {
		const idx_t offset(seq_idx[rid].offset);
		get_sequence(offset + from, offset + to, forward, seq);
	}
	static void decode_and_append_sequence(std::string& s, const char* const seq, idx_t i, const idx_t end_i) {
		s.reserve(s.size() + end_i - i);
		for (; i < end_i; ++i) {
			s += "ACGT"[static_cast<int>(seq[i])];
		}
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
	static void set_char(u1_t* const p, const idx_t idx, const u1_t c) {
		p[idx >> 2] |= c << ((~idx & 3) << 1);
	}
	static u1_t get_char(const u1_t* const p, const idx_t idx) {
		return p[idx >> 2] >> ((~idx & 3) << 1) & 3;
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
		return seq_idx[rid].offset;
	}
	idx_t seq_size(const idx_t rid) const {
		return seq_idx[rid].size;
	}
	void destroy() {
		if (pac) {
			safe_free(pac);
		}
		seq_idx.clear();
		db_size = max_db_size = 0;
	}
	void clear() {
		seq_idx.clear();
		db_size = 0;
		memset(pac, 0, (max_db_size + 3) / 4);
	}
	static void generate_pac_name(const char* const prefix, std::string& ret) {
		ret = prefix;
		ret += ".pac";
	}
	static void generate_idx_name(const char* const prefix, std::string& ret) {
		ret = prefix;
		ret += ".idx";
	}
	idx_t offset_to_rid(const idx_t offset) const;
	void add_one_seq(const Sequence& seq);
	void add_one_seq(const char* seq, const idx_t size);
	static void dump_pac(const u1_t* p, idx_t size, const char* path);
	static u1_t* load_pac(const char* path, idx_t& size);
	static void dump_idx(const PODArray<SeqIndex>& idx_list, const char* path);
	static void load_idx(const char* path, PODArray<SeqIndex>& idx_list);
	void dump_packed_db(const char* path) const;
	void load_packed_db(const char* path);
	static void pack_fasta_db(const char* fasta, const char* output_prefix, const idx_t min_size);
	void load_fasta_db(const char* fasta);
    private:
	u1_t* pac;
	idx_t db_size;
	idx_t max_db_size;
	PODArray<SeqIndex> seq_idx;
};

#endif // PACKED_DB_H
