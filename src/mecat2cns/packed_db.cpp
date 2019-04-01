#include "packed_db.h"

#include <fstream>	// ifstream, ofstream
#include <set>		// set<>
#include <stdio.h>	// rename()
#include <string>
#include <strings.h>	// bzero()
#include <sys/stat.h>	// stat(), struct stat
#include <unistd.h>	// unlink()

#include "../common/defs.h"
#include "../common/fasta_reader.h"
#include "../common/alignment.h"

void PackedDB::add_one_seq(const Sequence& seq) {
	SeqIndex si;
	si.file_offset = -1;
	si.memory_offset = db_size;
	si.size = seq.size();
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
	unsigned int rand_char(-1);	// spread out unknown sequence in a repeatable fashion
	for (idx_t i(0); i < si.size; ++i, ++db_size) {
		const u1_t c(table[static_cast<int>(org_seq[i])]);
		set_char(db_size, c < 4 ? c : ++rand_char & 3);
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

// XXX - make this restartable
void PackedDB::convert_fasta_to_db(const std::string& fasta, const std::string& output_prefix, const idx_t min_size) {
	DynamicTimer dtimer(__func__);
	u1_t buffer[MAX_SEQ_SIZE];
	const u1_t* const et(get_dna_encode_table());
	FastaReader fr(fasta.c_str());
	std::string filename(generate_pac_name(output_prefix));
	unlink(filename.c_str());
	filename += ".tmp";
	std::ofstream pout;
	open_fstream(pout, filename.c_str(), std::ios::out | std::ios::binary);
	std::streambuf* psb(pout.rdbuf());
	filename = generate_idx_name(output_prefix);
	unlink(filename.c_str());
	filename += ".tmp";
	std::ofstream iout;
	open_fstream(iout, filename.c_str(), std::ios::out);
	Sequence read;
	unsigned int rand_char(-1);	// spread out unknown sequence in a repeatable fashion
	off_t file_offset(0);
	for (;;) {
		const idx_t rsize(fr.read_one_seq(read));
		if (rsize == -1) {
			break;
		}
		if (rsize < min_size) {
			// can't skip entries in index or read ids won't
			// match ones from candidates, but put dummy size
			// in ones we don't use so we error out if they're
			// somehow used
			iout << file_offset << "\t0\n";
			continue;
		}
		Sequence::str_t& s(read.sequence());
		const idx_t rbytes((rsize + 3) / 4);
		// set_char uses | to set bits, so clear first
		bzero(buffer, rbytes);
		for (idx_t i(0); i < rsize; ++i) {
			const u1_t c(et[static_cast<int>(s[i])]);
			set_char(buffer, i, c < 4 ? c : ++rand_char & 3);
		}
		sb_write(psb, buffer, rbytes);
		iout << file_offset << "\t" << rsize << "\n";
		file_offset += rbytes;
	}
	close_fstream(pout);
	close_fstream(iout);
	filename = generate_pac_name(output_prefix);
	std::string tmp(filename + ".tmp");
	if (rename(tmp.c_str(), filename.c_str()) == -1) {
		ERROR("Could not rename tmp database file");
	}
	filename = generate_idx_name(output_prefix);
	tmp = filename + ".tmp";
	if (rename(tmp.c_str(), filename.c_str()) == -1) {
		ERROR("Could not rename tmp database index file");
	}
}

void PackedDB::open_db(const std::string& path, const idx_t size) {
	destroy();
	std::string filename(generate_pac_name(path));
	struct stat buf;
	if (stat(filename.c_str(), &buf) == -1) {
		ERROR("Could not stat fasta db file: %s", filename.c_str());
	}
	max_db_size = std::min(idx_t(buf.st_size), size);
	if (max_db_size > 0) {
		safe_calloc(pac, u1_t, max_db_size);
	}
	open_fstream(pstream, filename.c_str(), std::ios::in);
	filename = generate_idx_name(path);
	std::ifstream index;
	open_fstream(index, filename.c_str(), std::ios::in);
	SeqIndex si;
	if (max_db_size == buf.st_size && max_db_size) {	// read it all!
		while (index >> si.file_offset >> si.size) {
			si.memory_offset = si.file_offset * 4;
			seq_idx.push_back(si);
		}
		close_fstream(index);
		std::streambuf* sb(pstream.rdbuf());
		sb_read(sb, pac, max_db_size);
		close_fstream(pstream);
	} else {
		si.memory_offset = -1;
		while (index >> si.file_offset >> si.size) {
			seq_idx.push_back(si);
		}
		close_fstream(index);
	}
}

idx_t PackedDB::load_reads(const ExtensionCandidate* const ec_list, const idx_t nec) {
	if (!pstream.is_open()) {	// all in memory already
		return nec;
	}
	DynamicTimer dtimer(__func__);
	// not sure if it's faster to just read from disk and not
	// worry about obsolete memory_offset values, or if it would be
	// better to shift any reads in memory to bottom and then read in
	// new ones (which would require going through entire index and
	// would leave reads out of order)
	//
	// get set of reads to read in (most that will fit in memory);
	// use sorted set to speed reading (below)
	std::set<idx_t> read_ids;
	idx_t i(0), total_size(0);
	while (i != nec) {
		// find all alignments for a given read
		// (use set to handle duplicate qids)
		const idx_t sid(ec_list[i].sid);
		const idx_t start(i);
		std::set<idx_t> my_ids;
		my_ids.insert(sid);
		my_ids.insert(ec_list[i].qid);
		for (++i; i != nec && ec_list[i].sid == sid; ++i) {
			my_ids.insert(ec_list[i].qid);
		}
		// find size of new read additions
		idx_t size(0);
		std::set<idx_t>::const_iterator a(my_ids.begin());
		const std::set<idx_t>::const_iterator end_a(my_ids.end());
		for (; a != end_a; ++a) {
			if (read_ids.find(*a) == read_ids.end()) {
				size += (seq_idx[*a].size + 3) / 4;
			}
		}
		// if memory is limited, see if we've hit the limit
		if (max_db_size > 0 && total_size + size > max_db_size) {
			i = start;
			break;
		}
		total_size += size;
		read_ids.insert(my_ids.begin(), my_ids.end());
	}
	LOG(stderr, "using %ld bytes for %lu reads, %ld aligns (out of %ld)", total_size, read_ids.size(), i, nec);
	if (max_db_size == 0) {
		max_db_size = total_size;
		safe_calloc(pac, u1_t, max_db_size);
	}
	// now read in the reads
	std::set<idx_t>::const_iterator a(read_ids.begin());
	const std::set<idx_t>::const_iterator end_a(read_ids.end());
	idx_t pos(0);
	off_t offset(-1);
	for (; a != end_a; ++a) {
		SeqIndex& si(seq_idx[*a]);
		// don't seek if we're already in position
		if (offset != si.file_offset && !pstream.seekg(si.file_offset)) {
			ERROR("Error seeking on fasta db");
		}
		const idx_t bytes((si.size + 3) / 4);
		if (!pstream.read((char*)pac + pos, bytes)) {
			ERROR("Error reading fasta db");
		}
		si.memory_offset = pos * 4;
		pos += bytes;
		offset = si.file_offset + bytes;
	}
	if (i == nec) {
		pstream.close();
	}
	return i;
}

#if 0
// these routines aren't used anywhere
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
	const idx_t n(idx_list.size());
	for (idx_t i(0); i < n; ++i) {
		out << idx_list[i].memory_offset << "\t" << idx_list[i].size << "\n";
	}
	close_fstream(out);
}

void PackedDB::load_idx(const char* const path, PODArray<SeqIndex>& idx_list) {
	idx_list.clear();
	std::ifstream in;
	open_fstream(in, path, std::ios::in);
	SeqIndex si;
	si.file_offset = -1;
	while (in >> si.memory_offset >> si.size) {
		idx_list.push_back(si);
	}
	close_fstream(in);
}

void PackedDB::dump_packed_db(const char* const path) const {
	std::string(generate_pac_name(path));
	dump_pac(pac, db_size, n.c_str());
	n = generate_idx_name(path);
	dump_idx(seq_idx, n.c_str());
}

void PackedDB::load_packed_db(const char* const path) {
	std::string n(generate_pac_name(path));
	pac = load_pac(n.c_str(), db_size);
	max_db_size = db_size;
	n = generate_idx_name(path);
	load_idx(n.c_str(), seq_idx);
}

void PackedDB::pack_fasta_db(const char* const path, const char* const output_prefix, const idx_t min_size) {
	u1_t buffer[MAX_SEQ_SIZE];
	const u1_t* const et(get_dna_encode_table());
	FastaReader fr(path);
	std::string filename(generate_pac_name(output_prefix));
	std::ofstream pout;
	open_fstream(pout, filename.c_str(), std::ios::out | std::ios::binary);
	std::streambuf* psb(pout.rdbuf());
	filename = generate_idx_name(output_prefix);
	std::ofstream iout;
	open_fstream(iout, filename.c_str(), std::ios::out);
	Sequence read;
	unsigned int rand_char(-1);	// spread out unknown sequence in a repeatable fashion
	idx_t count(0), tsize(0);
	off_t file_offset(0);
	for (;;) {
		idx_t rsize(fr.read_one_seq(read));
		if (rsize == -1) {
			break;
		} else if (rsize < min_size) {
			continue;
		}
		Sequence::str_t& s(read.sequence());
		// set_char uses | to set bits, so clear first
		bzero(buffer, (rsize + 3) / 4);
		for (idx_t i(0); i < rsize; ++i) {
			const u1_t c(et[static_cast<int>(s[i])]);
			set_char(buffer, i, c < 4 ? c : ++rand_char & 3);
		}
		iout << file_offset << "\t" << rsize << "\n";
		tsize += rsize;
		rsize = (rsize + 3) / 4;
		sb_write(psb, buffer, rsize);
		file_offset += rsize * 4;
		++count;
	}
	sb_write(psb, &file_offset, sizeof(idx_t));
	close_fstream(pout);
	close_fstream(iout);
	LOG(stdout, "pack %ld reads, totally %ld residues", count, tsize);
}

void PackedDB::add_one_seq(const char* const seq, const idx_t size) {
	SeqIndex si;
	si.file_offset = -1;
	si.memory_offset = db_size;
	si.size = size;
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
	unsigned int rand_char(-1);	// spread out unknown sequence in a repeatable fashion
	for (idx_t i(0); i < si.size; ++i, ++db_size) {
		const u1_t c(table[static_cast<int>(seq[i])]);
		set_char(db_size, c < 4 ? c : ++rand_char & 3);
	}
}

// convert offset into read id
idx_t PackedDB::offset_to_rid(const idx_t offset) const {
	if (offset >= db_size) {
		return -1;
	}
	idx_t left(0), mid(0), right(seq_idx.size());
	while (left < right) {
		mid = (left + right) >> 1;
		if (offset < seq_idx[mid].memory_offset) {
			right = mid;
		} else if (mid == seq_idx.size() - 1) {
			break;
		} else if (offset < seq_idx[mid + 1].memory_offset) {
			break;
		} else {
			left = mid + 1;
		}
	}
	return mid;
}
#endif
