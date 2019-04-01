#include "packed_db.h"

#include <fstream>	// ifstream, ofstream
#include <set>		// set<>
#include <stdio.h>	// rename()
#include <string>
#include <strings.h>	// bzero()
#include <sys/stat.h>	// stat(), struct stat
#include <unistd.h>	// F_OK, unlink()

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

static int check_conversion_restart(const std::string& ckpt_file, off_t& fasta_offset, off_t& pac_offset, off_t& index_offset, unsigned int& rand_char) {
	std::ifstream in(ckpt_file.c_str());
	if (!in) {
		return 0;
	}
	in >> fasta_offset >> pac_offset >> index_offset >> rand_char;
	if (!in) {
		ERROR("Read error while restoring checkpoint from %s", ckpt_file.c_str());
	}
	return 1;
}

static void checkpoint_conversion(const std::string& ckpt_file, const std::string& ckpt_file_tmp, const off_t fasta_offset, const off_t pac_offset, const off_t index_offset, const unsigned int rand_char) {
	std::ofstream out(ckpt_file_tmp.c_str());
	if (!out) {
		LOG(stderr, "Checkpoint failed: couldn't open %s", ckpt_file_tmp.c_str());
		return;
	}
	out << fasta_offset << " " << pac_offset << " " << index_offset << " " << rand_char << "\n";
	if (!out) {
		LOG(stderr, "Checkpoint failed: write failed: %s", ckpt_file_tmp.c_str());
		return;
	}
	out.close();
	if (rename(ckpt_file_tmp.c_str(), ckpt_file.c_str()) == -1) {
		LOG(stderr, "Checkpoint failed: rename failed: %s", ckpt_file.c_str());
	}
}

void PackedDB::convert_fasta_to_db(const std::string& fasta, const std::string& output_prefix, const idx_t min_size) {
	const std::string pac_name(output_prefix + ".pac");
	const std::string index_name(output_prefix + ".idx");
	// see if we already did this
	if (access(pac_name.c_str(), F_OK) == 0 && access(index_name.c_str(), F_OK) == 0) {
		return;
	}
	DynamicTimer dtimer(__func__);
	u1_t buffer[MAX_SEQ_SIZE];
	const u1_t* const et(get_dna_encode_table());
	FastaReader fr(fasta.c_str());
	const std::string pac_name_tmp(pac_name + ".tmp");
	const std::string index_name_tmp(index_name + ".tmp");
	const std::string ckpt_name(pac_name + ".ckpt");
	const std::string ckpt_name_tmp(ckpt_name + ".tmp");
	std::ofstream pout, iout;
	unsigned int rand_char(-1);	// spread out unknown sequence in a repeatable fashion
	off_t file_offset(0), fasta_offset, index_offset;
	if (check_conversion_restart(ckpt_name, fasta_offset, file_offset, index_offset, rand_char)) {
		// don't truncate on restart
		open_fstream(pout, pac_name_tmp.c_str(), std::ios::out | std::ios::in | std::ios::binary);
		pout.seekp(file_offset);
		open_fstream(iout, index_name_tmp.c_str(), std::ios::out | std::ios::in);
		iout.seekp(index_offset);
		fr.seekg(fasta_offset);
	} else {
		open_fstream(pout, pac_name_tmp.c_str(), std::ios::out | std::ios::binary);
		open_fstream(iout, index_name_tmp.c_str(), std::ios::out);
	}
	Sequence read;
	time_t next_checkpoint_time(time(0) + 300);
	for (;;) {
		const idx_t rsize(fr.read_one_seq(read));
		if (rsize == -1) {
			break;
		}
		if (rsize < min_size) {
			// can't skip entries in index or read ids won't
			// match ones from candidates, so put dummy size
			// in ones we don't use
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
		if (!pout.write((char*)buffer, rbytes)) {
			ERROR("Write error to file %s", pac_name_tmp.c_str());
		}
		iout << file_offset << "\t" << rsize << "\n";
		file_offset += rbytes;
		if (time(0) >= next_checkpoint_time) {
			checkpoint_conversion(ckpt_name, ckpt_name_tmp, fr.tellg(), file_offset, iout.tellp(), rand_char);
			next_checkpoint_time = time(0) + 300;
		}
	}
	checkpoint_conversion(ckpt_name, ckpt_name_tmp, fr.tellg(), file_offset, iout.tellp(), rand_char);
	close_fstream(pout);
	close_fstream(iout);
	if (rename(pac_name_tmp.c_str(), pac_name.c_str()) == -1) {
		ERROR("Could not rename tmp database file");
	}
	if (rename(index_name_tmp.c_str(), index_name.c_str()) == -1) {
		ERROR("Could not rename tmp database index file");
	}
	unlink(ckpt_name.c_str());
}

void PackedDB::open_db(const std::string& path, const idx_t size) {
	destroy();
	const std::string pac_name(path + ".pac");
	struct stat buf;
	if (stat(pac_name.c_str(), &buf) == -1) {
		ERROR("Could not stat fasta db file: %s", pac_name.c_str());
	}
	max_db_size = std::min(idx_t(buf.st_size), size);
	if (max_db_size) {
		safe_calloc(pac, u1_t, max_db_size);
	}
	open_fstream(pstream, pac_name.c_str(), std::ios::in);
	const std::string index_name(path + ".idx");
	std::ifstream index;
	open_fstream(index, index_name.c_str(), std::ios::in);
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
	// would leave reads out of order); past experience says the former
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
		if (max_db_size && total_size + size > max_db_size) {
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
