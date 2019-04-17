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
	unsigned int rand_char(0);	// spread out unknown sequence in a repeatable fashion
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

static int check_conversion_restart(const std::string& ckpt_file, off_t& fasta_offset, off_t& pac_offset, off_t& index_offset, unsigned int& rand_char, size_t& read_count) {
	std::ifstream in(ckpt_file.c_str());
	if (!in) {
		return 0;
	}
	in >> fasta_offset >> pac_offset >> index_offset >> rand_char >> read_count;
	if (!in) {
		ERROR("Read error while restoring checkpoint from %s", ckpt_file.c_str());
	}
	return 1;
}

static void checkpoint_conversion(const std::string& ckpt_file, const std::string& ckpt_file_tmp, const off_t fasta_offset, const off_t pac_offset, const off_t index_offset, const unsigned int rand_char, size_t read_count) {
	std::ofstream out(ckpt_file_tmp.c_str());
	if (!out) {
		LOG(stderr, "Checkpoint failed: couldn't open %s", ckpt_file_tmp.c_str());
		return;
	}
	out << fasta_offset << " " << pac_offset << " " << index_offset << " " << rand_char << " " << read_count << "\n";
	if (!out) {
		LOG(stderr, "Checkpoint failed: write failed: %s", ckpt_file_tmp.c_str());
		return;
	}
	out.close();
	if (rename(ckpt_file_tmp.c_str(), ckpt_file.c_str()) == -1) {
		LOG(stderr, "Checkpoint failed: rename failed: %s", ckpt_file.c_str());
	}
}

size_t PackedDB::convert_fasta_to_db(const std::string& fasta, const std::string& output_prefix, const idx_t min_size) {
	const std::string pac_name(output_prefix + ".pac");
	const std::string index_name(output_prefix + ".idx");
	size_t read_count(0);
	// see if we already did this
	if (access(pac_name.c_str(), F_OK) == 0 && access(index_name.c_str(), F_OK) == 0) {
		// get number of reads from end of database
		std::ifstream pin;
		open_fstream(pin, pac_name.c_str(), std::ios::in);
		if (!pin.seekg(-sizeof(size_t), std::ios_base::end)) {
			ERROR("Could not seek to end of fasta db to get size\n");
		}
		if (!pin.read((char*)&read_count, sizeof(size_t))) {
			ERROR("Could not read fasta db to get size\n");
		}
		close_fstream(pin);
		return read_count;
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
	unsigned int rand_char(0);	// spread out unknown sequence in a repeatable fashion
	off_t pac_offset(0), fasta_offset, index_offset;
	if (check_conversion_restart(ckpt_name, fasta_offset, pac_offset, index_offset, rand_char, read_count)) {
		// don't truncate on restart
		open_fstream(pout, pac_name_tmp.c_str(), std::ios::out | std::ios::in | std::ios::binary);
		pout.seekp(pac_offset);
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
		++read_count;
		if (rsize < min_size) {
			// can't skip entries in index or read ids won't
			// match ones from candidates, so put dummy size
			// in ones we don't use
			iout << pac_offset << "\t0\n";
			if (!iout) {
				ERROR("Error writing to pac index file: %s", index_name_tmp.c_str());
			}
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
		iout << pac_offset << "\t" << rsize << "\n";
		if (!iout) {
			ERROR("Error writing to pac index file: %s", index_name_tmp.c_str());
		}
		pac_offset += rbytes;
		if (time(0) >= next_checkpoint_time) {
			checkpoint_conversion(ckpt_name, ckpt_name_tmp, fr.tellg(), pac_offset, iout.tellp(), rand_char, read_count);
			next_checkpoint_time = time(0) + 300;
		}
	}
	checkpoint_conversion(ckpt_name, ckpt_name_tmp, fr.tellg(), pac_offset, iout.tellp(), rand_char, read_count);
	if (!pout.write((char*)&read_count, sizeof(size_t))) {
		ERROR("Write error to file %s", pac_name_tmp.c_str());
	}
	close_fstream(pout);
	close_fstream(iout);
	if (rename(pac_name_tmp.c_str(), pac_name.c_str()) == -1) {
		ERROR("Could not rename tmp database file");
	}
	if (rename(index_name_tmp.c_str(), index_name.c_str()) == -1) {
		ERROR("Could not rename tmp database index file");
	}
	unlink(ckpt_name.c_str());
	return read_count;
}

void PackedDB::open_db(const std::string& path, const idx_t size) {
	destroy();
	const std::string pac_name(path + ".pac");
	open_fstream(pstream, pac_name.c_str(), std::ios::in);
	// get number of reads from end of database
	if (!pstream.seekg(-sizeof(size_t), std::ios_base::end)) {
		ERROR("Could not seek to end of fasta db to get size\n");
	}
	const idx_t file_size(pstream.tellg());
	max_db_size = size ? std::min(file_size, size) : file_size;
	if (max_db_size) {
		safe_calloc(pac, u1_t, max_db_size);
	}
	size_t read_count;
	if (!pstream.read((char*)&read_count, sizeof(size_t))) {
		ERROR("Could not read fasta db to get size\n");
	}
	seq_idx.reserve(read_count);
	pstream.clear();					// clear eof
	if (!pstream.seekg(0, std::ios_base::beg)) {
		ERROR("Could not reset to start of fast db file\n");
	}
	const std::string index_name(path + ".idx");
	std::ifstream index;
	open_fstream(index, index_name.c_str(), std::ios::in);
	SeqIndex si;
	if (max_db_size == file_size && max_db_size) {		// read it all!
		while (index >> si.file_offset >> si.size) {
			si.memory_offset = si.file_offset * 4;
			seq_idx.push_back(si);
		}
		close_fstream(index);
		if (!pstream.read((char*)pac, max_db_size)) {
			ERROR("Error reading fasta database\n");
		}
		close_fstream(pstream);
	} else {
		si.memory_offset = -1;
		while (index >> si.file_offset >> si.size) {
			seq_idx.push_back(si);
		}
		close_fstream(index);
	}
}

const char* PackedDB::load_read(const idx_t read_id) {
	const SeqIndex& si(seq_idx[read_id]);
	if (!pstream.is_open()) {	// all in memory already
		return (char*)pac + si.memory_offset / 4;
	}
	if (!pstream.seekg(si.file_offset)) {
		ERROR("Error seeking on fasta db");
	}
	const idx_t bytes((si.size + 3) / 4);
	if (!pstream.read((char*)pac, bytes)) {
		ERROR("Error reading fasta db");
	}
	return (char*)pac;
}

idx_t PackedDB::load_reads(const ExtensionCandidateCompressed* const ec_list, const idx_t nec) {
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

void PackedDB::create_index(const std::string& output_prefix, const std::vector<std::pair<idx_t, idx_t> >& index) {
	const std::string index_name(output_prefix + ".idx");
	const std::string index_name_tmp(index_name + ".tmp");
	std::ofstream out;
	open_fstream(out, index_name_tmp.c_str(), std::ios::out);
	std::vector<std::pair<idx_t, idx_t> >::const_iterator a(index.begin());
	const std::vector<std::pair<idx_t, idx_t> >::const_iterator end_a(index.end());
	for (; a != end_a; ++a) {
		out << a->first << "\t" << a->second << "\n";
		if (!out) {
			ERROR("Error writing to pac index file: %s", index_name_tmp.c_str());
		}
	}
	close_fstream(out);
	if (rename(index_name_tmp.c_str(), index_name.c_str()) == -1) {
		ERROR("Error renaming pac index file: %s", index_name_tmp.c_str());
	}
}

void PackedDB::read_index(const std::string& output_prefix, std::vector<std::pair<idx_t, idx_t> >& index) {
	const std::string pac_name(output_prefix + ".pac");
	std::ifstream in(pac_name.c_str());
	// pre-allocate index if possible
	if (in.is_open()) {
		// get number of reads from end of database
		if (!in.seekg(-sizeof(size_t), std::ios_base::end)) {
			ERROR("Could not seek to end of fasta db to get size\n");
		}
		size_t read_count;
		if (!in.read((char*)&read_count, sizeof(size_t))) {
			ERROR("Could not read fasta db to get size\n");
		}
		in.close();
		index.reserve(read_count);
	}
	const std::string index_name(output_prefix + ".idx");
	open_fstream(in, index_name.c_str(), std::ios::in);
	idx_t i, j;
	while (in >> i >> j) {
		index.push_back(std::make_pair(i, j));
	}
}

void PackedDB::read_sizes(const std::string& output_prefix, std::vector<idx_t>& sizes) {
	const std::string pac_name(output_prefix + ".pac");
	std::ifstream in(pac_name.c_str());
	// pre-allocate index if possible
	if (in.is_open()) {
		// get number of reads from end of database
		if (!in.seekg(-sizeof(size_t), std::ios_base::end)) {
			ERROR("Could not seek to end of fasta db to get size\n");
		}
		size_t read_count;
		if (!in.read((char*)&read_count, sizeof(size_t))) {
			ERROR("Could not read fasta db to get size\n");
		}
		in.close();
		sizes.reserve(read_count);
	}
	const std::string index_name(output_prefix + ".idx");
	open_fstream(in, index_name.c_str(), std::ios::in);
	idx_t i, j;
	while (in >> i >> j) {
		sizes.push_back(j);
	}
}
