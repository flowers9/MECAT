#ifndef _READS_CORRECTION_AUX_H
#define _READS_CORRECTION_AUX_H

#include <vector>
#include <cstring>
#include <unistd.h>

#include "dw.h"
#include "packed_db.h"
#include "options.h"

struct CnsTableItem {
	char base;
	uint1 mat_cnt;
	uint1 ins_cnt;
	uint1 del_cnt;
	CnsTableItem() : base('N'), mat_cnt(0), ins_cnt(0), del_cnt(0) { }
};

#define MAX_CNS_OVLPS 100

class MappingRange {
    public:
	int start, end;
	explicit MappingRange() : start(0), end(0) { }
	explicit MappingRange(const int s, const int e) : start(s), end(e) { }
	~MappingRange() { }
};

class CnsAln : public MappingRange {
    public:
	explicit CnsAln() : aln_idx(0) { }
	explicit CnsAln(const int i, const int j, const int k, const std::string& q, const std::string& s) : MappingRange(i, j), aln_idx(k), qaln(s), saln(q) { }
	~CnsAln() { }
	// don't know why this skips the first basepair
	int retrieve_aln_subseqs(const int sb, const int se, std::string& qstr, std::string& tstr, int& sb_out) {
		const int aln_size(static_cast<int>(saln.size()) - 1);
		if (se <= start || sb >= end || aln_idx >= aln_size) {
			return 0;
		}
		sb_out = std::max(start, sb);
		qstr.clear();
		tstr.clear();
		while (start < sb && aln_idx < aln_size) {
			++aln_idx;
			if (saln[aln_idx] != GAP) {
				++start;
			}
		}
		qstr += qaln[aln_idx];
		tstr += saln[aln_idx];
		while (start < se && aln_idx < aln_size) {
			++aln_idx;
			if (saln[aln_idx] != GAP) {
				++start;
			}
			qstr += qaln[aln_idx];
			tstr += saln[aln_idx];
		}
		return 1;
	}
    private:
	int aln_idx;
	std::string qaln, saln;
};

class CnsAlns {
    public:
	CnsAlns() { }
	~CnsAlns() { }
	void clear() {
		cns_alns_.clear();
	}
	size_t num_alns() const {
		return cns_alns_.size();
	}
	std::vector<CnsAln>::iterator begin() {
		return cns_alns_.begin();
	}
	std::vector<CnsAln>::const_iterator end() const {
		return cns_alns_.end();
	}
	void add_aln(const int soff, const int send, const std::string& qstr, const std::string& tstr) {
		r_assert(qstr.size() == tstr.size());
		cns_alns_.push_back(CnsAln(soff, send, 0, qstr, tstr));
	}
	void get_mapping_ranges(std::vector<MappingRange>& ranges) const {
		ranges.clear();
		ranges.reserve(cns_alns_.size());
		for (size_t i(0); i < cns_alns_.size(); ++i) {
			ranges.push_back(cns_alns_[i]);
		}
	}
    private:
	std::vector<CnsAln> cns_alns_;
};

// 1k seems to work a bit better than 10k - perhaps less time waiting for
// another thread to finish writing?
#define MAX_CNS_RESULTS 1000

struct CmpExtensionCandidateBySidAndScore {
	bool operator()(const ExtensionCandidate& a, const ExtensionCandidate& b) {
		if (a.sid != b.sid) {			// primary sort
			return a.sid < b.sid;
		} else if (a.score != b.score) {	// secondary sort
			return b.score < a.score;	// process best ones first
		} else if (a.qid != b.qid) {
			return a.qid < b.qid;
		} else if (a.qext != b.qext) {
			return a.qext < b.qext;
		} else if (a.sext != b.sext) {		// tertiary sort
			return a.sext < b.sext;		// make sorting consistent
		} else if (a.soff != b.soff) {
			return a.soff < b.soff;
		} else if (a.qoff != b.qoff) {
			return a.qoff < b.qoff;
		} else if (a.send != b.send) {
			return a.send < b.send;
		} else if (a.qend != b.qend) {
			return a.qend < b.qend;
		} else {
			return a.qdir < b.qdir;
		}
		// sdir, qsize, ssize are all the same by this point
	}
};

struct CmpExtensionCandidateCompressedBySidAndScore {
	bool operator()(const ExtensionCandidateCompressed& a, const ExtensionCandidateCompressed& b) {
		if (a.sid != b.sid) {			// primary sort
			return a.sid < b.sid;		// for splitting up in allocate_ecs()
		} else if (a.score != b.score) {	// secondary sort
			return b.score < a.score;	// process best ones first
		} else if (a.qid != b.qid) {
			return a.qid < b.qid;
		} else if (a.qext() != b.qext()) {
			return a.qext() < b.qext();
		} else if (a.sext != b.sext) {		// tertiary sort
			return a.sext < b.sext;		// make sorting consistent
		} else {
			return a.qdir() < b.qdir();
		}
	}
};

class ConsensusPerThreadData {
    public:
	// this is ExtensionCandidate for m4 runs, ExtensionCandidateCompressed
	// for candidate runs (to reduce memory usage)
	void* candidates;
	// num_candidates, candidates initialized by allocate_ecs()
	// next_candidate initialized by ConsensusThreadData::restart()
	idx_t num_candidates, next_candidate;
	DiffRunningData drd;
	M5Record m5;
	CnsAlns cns_alns;
	std::vector<CnsTableItem> cns_table;
	std::vector<uint1> id_list;
	std::vector<CnsResult> cns_results;
	std::string query, target, qaln, saln;
    public:
	ConsensusPerThreadData() : drd(DiffRunningData(get_sw_parameters_small())) {
		// we'll definitely be seeing at least this much use,
		// so might as well preallocate
		cns_results.reserve(MAX_CNS_RESULTS);
	}
	~ConsensusPerThreadData() { }
};

class ConsensusThreadData {
    public:
	ReadsCorrectionOptions& rco;
	PackedDB& reads;
	std::ostream& out;
	ConsensusPerThreadData* data;
	pthread_mutex_t out_lock;
	idx_t ec_offset;
	// this doesn't work as a vector - all the pointers end up pointing
	// to the same values, and eventually it seg faults (possibly a
	// compiler optimization bug)
    public:
	ConsensusThreadData(ReadsCorrectionOptions& prco, PackedDB& r, std::ostream& output, const char* const input_file_name) : rco(prco), reads(r), out(output), data(new ConsensusPerThreadData[prco.num_threads]), ec_offset(0), last_thread_id_(-1), num_threads_written_(0) {
		done_file_ = input_file_name;
		done_file_ += ".done";
		ckpt_file_ = input_file_name;
		ckpt_file_ += ".ckpt";
		ckpt_file_tmp_ = ckpt_file_ + ".tmp";
		pthread_mutex_init(&out_lock, NULL);
		pthread_mutex_init(&id_lock_, NULL);
	}
	~ConsensusThreadData() {
		delete[] data;
		pthread_mutex_destroy(&out_lock);
		pthread_mutex_destroy(&id_lock_);
	}
	int get_thread_id() {
		pthread_mutex_lock(&id_lock_);
		const int tid(++last_thread_id_);
		pthread_mutex_unlock(&id_lock_);
		return tid;
	}
	void reset_threads() {
		last_thread_id_ = -1;
		for (int i(0); i < rco.num_threads; ++i) {
			data[i].next_candidate = 0;
		}
	}
	void write_buffer(const int tid, const idx_t i) {
		ConsensusPerThreadData& pdata(data[tid]);
		std::vector<CnsResult>::const_iterator a(pdata.cns_results.begin());
		const std::vector<CnsResult>::const_iterator end_a(pdata.cns_results.end());
		pthread_mutex_lock(&out_lock);
		for (; a != end_a; ++a) {
			out << ">" << a->id << "_" << a->range[0] << "_" << a->range[1] << "_" << a->seq.size() << "\n" << a->seq << "\n";
			if (!out) {
				ERROR("Error writing output");
			}
		}
		pdata.next_candidate = i;
		if (++num_threads_written_ >= rco.num_threads) {
			checkpoint();
			num_threads_written_ = 0;
		}
		pthread_mutex_unlock(&out_lock);
		pdata.cns_results.clear();
	}
	int restart(off_t& output_pos) {
		if (access(ckpt_file_.c_str(), F_OK) != 0) {
			for (int i(0); i < rco.num_threads; ++i) {
				data[i].next_candidate = 0;
			}
			return 0;
		}
		std::ifstream ckpt_in(ckpt_file_.c_str());
		if (!ckpt_in) {
			ERROR("Restart failed: could not open checkpoint file: %s", ckpt_file_.c_str());
		}
		ckpt_in >> rco.job_index >> output_pos >> ec_offset;
		if (!ckpt_in) {
			ERROR("Restart failed: could not read checkpoint file: %s", ckpt_file_.c_str());
		}
		for (int i(0); i < rco.num_threads; ++i) {
			ckpt_in >> data[i].next_candidate;
			if (!ckpt_in) {
				ERROR("Restart failed: could not read checkpoint file: %s", ckpt_file_.c_str());
			}
		}
		return 1;
	}
    private:
	void checkpoint() {
		std::ofstream ckpt_out(ckpt_file_tmp_.c_str());
		if (!ckpt_out) {
			LOG(stderr, "Checkpoint failed: couldn't open %s", ckpt_file_tmp_.c_str());
			return;
		}
		out.flush();
		ckpt_out << rco.job_index << " " << off_t(out.tellp()) << " " << ec_offset << "\n";
		if (!ckpt_out) {
			LOG(stderr, "Checkpoint failed: write failed: %s", ckpt_file_tmp_.c_str());
			return;
		}
		for (int i(0); i < rco.num_threads; ++i) {
			ckpt_out << data[i].next_candidate << "\n";
			if (!ckpt_out) {
				LOG(stderr, "Checkpoint failed: write failed: %s", ckpt_file_tmp_.c_str());
				return;
			}
		}
		ckpt_out.close();
		if (rename(ckpt_file_tmp_.c_str(), ckpt_file_.c_str()) == -1) {
			LOG(stderr, "Checkpoint failed: rename failed: %s", ckpt_file_.c_str());
		}
	}
    private:
	pthread_mutex_t id_lock_;
	int last_thread_id_, num_threads_written_;
	std::string done_file_, ckpt_file_, ckpt_file_tmp_;
};

void normalize_gaps(const std::string& qstr, const std::string& tstr, std::string& qnorm, std::string& tnorm, int push);

void allocate_ecs(ConsensusThreadData &data, ExtensionCandidate* ec_list, idx_t nec);
void allocate_ecs(ConsensusThreadData &data, ExtensionCandidateCompressed* ec_list, idx_t nec);

#endif // _READS_CORRECTION_AUX_H
