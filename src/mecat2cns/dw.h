#ifndef DW_H
#define DW_H

#include "../common/defs.h"	// idx_t
#include <math.h>	// ceil()
#include <string>	// string
#include <vector>	// vector<>

class Alignment {
    public:
	// size tracks actual buffer use
	int aln_q_e, aln_t_e, size;
	std::vector<char> q_aln_str, t_aln_str;
    public:
	explicit Alignment() { }
	~Alignment() { }
	void resize(const size_t max_size) {
		q_aln_str.resize(max_size);
		t_aln_str.resize(max_size);
	}
	void reset() {
		size = 0;
	}
};

class OutputStore {
    public:
	// these track actual buffer use - buffer use starts at
	// buffer_start, with left going down, and right going up
	int buffer_start, left_size, right_size;
	int query_start, query_end;
	int target_start, target_end;
	std::vector<char> q_buffer, t_buffer;
    public:
	explicit OutputStore() { }
	~OutputStore() { }
	// can't figure out how to pass this in on initialization
	void resize(const size_t max_size) {
		q_buffer.resize(max_size);
		t_buffer.resize(max_size);
	}
	void reset_buffer(const int i) {
		buffer_start = i;
		left_size = right_size = 0;
	}
};

struct DPathData {
	int x1, y1, x2, y2, pre_k;
	explicit DPathData() { }
	explicit DPathData(const int i, const int j, const int k, const int l, const int m) : x1(i), y1(j), x2(k), y2(l), pre_k(m) { }
	void set(const int i, const int j, const int k, const int l, const int m) {
		x1 = i;
		y1 = j;
		x2 = k;
		y2 = l;
		pre_k = m;
	}
};

struct DPathIndex {
	int d_offset, min_k;
	explicit DPathIndex() { }
	void set(const int i, const int j) {
		d_offset = i;
		min_k = j;
	}
};

struct PathPoint {
	int x, y;
	explicit PathPoint() { }
	void set(const int i, const int j) {
		x = i;
		y = j;
	}
};

// if the end of query/target is within this or less, extend the whole distance
#define SEGMENT_BORDER 100

// use static sizing here rather than dynamic, as dynamic is 50% slower in practice;
// that said, reserving the space instead of initializing it may be viable

class DiffRunningData {
    public:
	static const int segment_size = 500;		// 500 is "small", 1000 is "large"
	Alignment align;
	OutputStore result;
	std::vector<int> DynQ, DynT;
	std::vector<DPathData> d_path;
	std::vector<DPathIndex> d_path_index;
	std::vector<PathPoint> aln_path;
    public:
	explicit DiffRunningData() { }
	~DiffRunningData() { }
	// can't figure out how to pass these in on initialization
	void set_size(const double error_rate, const idx_t max_read_size) {
		const size_t max_extend_size(segment_size + SEGMENT_BORDER);
		// in Align(), k_offset = extend_size * 4 * error_rate,
		const size_t max_k_offset(ceil(max_extend_size * 4 * error_rate));
		// allocate largest first (approximately)
		d_path.resize(max_k_offset * (max_k_offset + 1) / 2);
		result.resize(max_read_size * 2);
		aln_path.resize(max_k_offset * 4);
		align.resize(max_extend_size * 2);
		DynQ.resize(max_k_offset * 2);
		DynT.resize(max_k_offset * 2);
		d_path_index.resize(max_k_offset);
	}
};

class M5Record {
    public:
	idx_t qoff, qend, soff, send;
	std::string qaln, saln;
	M5Record() { }
	~M5Record() { }
};

int GetAlignment(const std::string& query, int query_start, const std::string& target, int target_start, DiffRunningData& drd, M5Record& m5, double error_rate, int min_aln_size);

#endif  // DW_H
