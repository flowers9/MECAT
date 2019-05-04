#ifndef DW_H
#define DW_H

#include <string>	// string
#include <vector>	// vector<>
#include "../common/defs.h"	// idx_t

struct SW_Parameters {
	int segment_size;	// probably best if it's a multiple of 10
	int row_size, column_size, d_path_size;
	SW_Parameters(const int i, const int j, const int k, const int l) : segment_size(i), row_size(j), column_size(k), d_path_size(l) { }
};

inline SW_Parameters get_sw_parameters_small() {
	// 1000 instead of 500 for "large"
	return SW_Parameters(500, 4096, 4096, 500000);
}

class Alignment {
    public:
	// size tracks actual buffer use
	int aln_q_e, aln_t_e, size;
	// these are buffers we only expand
	std::vector<char> q_aln_str, t_aln_str;
    public:
	explicit Alignment() { }
	explicit Alignment(const size_t new_max_size) {
		q_aln_str.resize(new_max_size);
		t_aln_str.resize(new_max_size);
	}
	~Alignment() { }
	void reset() {
		size = 0;
	}
};

class OutputStore {
    public:
	// these track actual buffer use
	int buffer_start, left_size, right_size;
	int query_start, query_end;
	int target_start, target_end;
	// these are buffers that we only expand
	std::vector<char> q_buffer, t_buffer;
	// for the record, inserts are q_buffer == 4, deletes are t_buffer == 4;
	// matches/mismatches are pretty obvious
    public:
	explicit OutputStore() { }
	explicit OutputStore(const size_t new_max_size) {
		q_buffer.resize(new_max_size);
		t_buffer.resize(new_max_size);
	}
	~OutputStore() { }
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

class DiffRunningData {
    public:
	const int segment_size;
	Alignment align;
	OutputStore result;
	std::vector<int> DynQ, DynT;
	std::vector<DPathData> d_path;
	std::vector<DPathIndex> d_path_index;
	std::vector<PathPoint> aln_path;
    public:
	// work on more accurate sizing of these buffers
	explicit DiffRunningData(const SW_Parameters& swp) : segment_size(swp.segment_size), align((swp.segment_size + 100) * 2), result(MAX_SEQ_SIZE * 2), DynQ(swp.row_size), DynT(swp.column_size), d_path(swp.d_path_size), d_path_index(swp.segment_size * 2), aln_path(swp.segment_size * 4) { }
	~DiffRunningData() { }
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
