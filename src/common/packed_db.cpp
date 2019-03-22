#include "packed_db.h"

#include <fstream>
#include <string>

#include "defs.h"
#include "fasta_reader.h"

void PackedDB::dump_pac(const u1_t* const p, const idx_t size, const char* const path) {
	std::ofstream out;
	open_fstream(out, path, std::ios::out | std::ios::binary);
	std::streambuf* sb(out.rdbuf());
	sb_write(sb, p, (size + 3) / 4);
	sb_write(sb, &size, sizeof(idx_t));
	close_fstream(out);
}

u1_t* PackedDB::load_pac(const char* path, idx_t& size) {
	std::ifstream in;
	open_fstream(in, path, std::ios::in | std::ios::binary);
	std::streambuf* sb(in.rdbuf());
	in.seekg(-sizeof(idx_t), std::ios::end);
	sb_read(sb, &size, sizeof(idx_t));
	in.seekg(0, std::ios::beg);
	u1_t* p;
	safe_calloc(p, u1_t, (size + 3) / 4);
	sb_read(sb, p, (size + 3) / 4);
	close_fstream(in);
	return p;
}

void PackedDB::dump_idx(const PODArray<SeqIndex>& idx_list, const char* const path) {
	std::ofstream out;
	open_fstream(out, path, std::ios::out);
	idx_t i(0), n(idx_list.size());
	for (i = 0; i < n; ++i) {
		out << idx_list[i].offset << "\t" << idx_list[i].size << "\n";
	}
	close_fstream(out);
}

void PackedDB::load_idx(const char* const path, PODArray<SeqIndex>& idx_list) {
	idx_list.clear();
	std::ifstream in;
	open_fstream(in, path, std::ios::in);
	SeqIndex si;
	while (in >> si.offset >> si.size) {
		idx_list.push_back(si);
	}
	close_fstream(in);
}

void PackedDB::dump_packed_db(const char* const path) const {
	std::string n;
	generate_pac_name(path, n);
	dump_pac(pac, db_size, n.c_str());
	generate_idx_name(path, n);
	dump_idx(seq_idx, n.c_str());
}

void PackedDB::load_packed_db(const char* const path) {
	std::string n;
	generate_pac_name(path, n);
	pac = load_pac(n.c_str(), db_size);
	max_db_size = db_size;
	generate_idx_name(path, n);
	load_idx(n.c_str(), seq_idx);
}

void PackedDB::pack_fasta_db(const char* const path, const char* const output_prefix, const idx_t min_size) {
	u1_t* buffer;
	safe_malloc(buffer, u1_t, MAX_SEQ_SIZE);
	const u1_t* const et(get_dna_encode_table());
	FastaReader fr(path);
	std::string n;
	generate_pac_name(output_prefix, n);
	std::ofstream pout;
	open_fstream(pout, n.c_str(), std::ios::out | std::ios::binary);
	std::streambuf* psb(pout.rdbuf());
	generate_idx_name(output_prefix, n);
	std::ofstream iout;
	open_fstream(iout, n.c_str(), std::ios::out);
	Sequence read;
	idx_t count(0), tsize(0);
	for (;;) {
		idx_t rsize(fr.read_one_seq(read));
		if (rsize == -1) {
			break;
		} else if (rsize < min_size) {
			continue;
		}
		Sequence::str_t& s(read.sequence());
		memset(buffer, 0, MAX_SEQ_SIZE);
		for (idx_t i(0); i < rsize; ++i) {
			const u1_t c(et[static_cast<int>(s[i])]);
			set_char(buffer, i, c < 3 ? c : 3);
		}
		iout << tsize << "\t" << rsize << "\n";
		rsize = (rsize + 3) / 4;
		sb_write(psb, buffer, rsize);
		tsize += rsize * 4;
		++count;
	}
	sb_write(psb, &tsize, sizeof(idx_t));
	close_fstream(pout);
	close_fstream(iout);
	safe_free(buffer);
	LOG(stdout, "pack %lld reads, totally %lld residues", (long long)count, (long long)tsize);
}

void PackedDB::add_one_seq(const Sequence& seq) {
	SeqIndex si;
	si.size = seq.size();
	si.offset = db_size;
	seq_idx.push_back(si);
	const idx_t needed_size(db_size + si.size);
	if (max_db_size < needed_size) {
		idx_t new_size(max_db_size > 1024 ? max_db_size : 1024);
		for (; new_size < needed_size; new_size *= 2) { }
		u1_t* new_pac(NULL);
		safe_calloc(new_pac, u1_t, (new_size + 3) / 4);
		memcpy(new_pac, pac, (db_size + 3) / 4);
		safe_free(pac);
		pac = new_pac;
		max_db_size = new_size;
	}
	const Sequence::str_t& org_seq(seq.sequence());
	const u1_t* const table(get_dna_encode_table());
	for (idx_t i(0); i < si.size; ++i, ++db_size) {
		const u1_t c(table[static_cast<int>(org_seq[i])]);
		set_char(db_size, c < 3 ? c : 3);
	}
}

void PackedDB::add_one_seq(const char* const seq, const idx_t size) {
	SeqIndex si;
	si.size = size;
	si.offset = db_size;
	seq_idx.push_back(si);
	const idx_t needed_size(db_size + si.size);
	if (max_db_size < needed_size) {
		idx_t new_size(max_db_size > 1024 ? max_db_size : 1024);
		for (; new_size < needed_size; new_size *= 2) { }
		u1_t* new_pac(NULL);
		safe_calloc(new_pac, u1_t, (new_size + 3) / 4);
		memcpy(new_pac, pac, (db_size + 3) / 4);
		safe_free(pac);
		pac = new_pac;
		max_db_size = new_size;
	}
	const u1_t* const table(get_dna_encode_table());
	for (idx_t i(0); i < si.size; ++i, ++db_size) {
		const u1_t c(table[static_cast<int>(seq[i])]);
		set_char(db_size, c < 3 ? c : 3);
	}
}

void PackedDB::load_fasta_db(const char* const dbname) {
	DynamicTimer dtimer(__func__);
	FastaReader freader(dbname);
	Sequence seq;
	while (freader.read_one_seq(seq) != -1) {
		add_one_seq(seq);
	}
}

idx_t PackedDB::offset_to_rid(const idx_t offset) const {
	if (offset >= db_size) {
		return -1;
	}
	idx_t left(0), mid(0), right(seq_idx.size());
	while (left < right) {
		mid = (left + right) >> 1;
		if (offset < seq_idx[mid].offset) {
			right = mid;
		} else if (mid == seq_idx.size() - 1) {
			break;
		} else if (offset < seq_idx[mid + 1].offset) {
			break;
		} else {
			left = mid + 1;
		}
	}
	return mid;
}
