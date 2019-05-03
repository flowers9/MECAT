#ifndef DW_H
#define DW_H

#include <string>	// string
#include <vector>	// vector<>
#include "../common/defs.h"	// idx_t

struct SW_Parameters {
	int segment_size;	// probably best if it's a multiple of 10
	int row_size, column_size;
	SW_Parameters(const int i, const int j, const int k) : segment_size(i), row_size(j), column_size(k) { }
};

inline SW_Parameters get_sw_parameters_small() {
	// 1000 instead of 500 for "large"
	return SW_Parameters(500, 4096, 4096);
}

class Alignment {
    public:
	// size tracks actual buffer use
	int aln_q_e, aln_t_e, size;
	// these are buffers we only expand
	std::vector<char> q_aln_str, t_aln_str;
    public:
	explicit Alignment() { }
	~Alignment() { }
	void reset(const size_t new_max_size) {
		size = 0;
		if (q_aln_str.size() < new_max_size) {
			q_aln_str.resize(new_max_size);
			t_aln_str.resize(new_max_size);
		}
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
	~OutputStore() { }
	void reset_buffer(const int i, const size_t new_max_size) {
		buffer_start = i;
		left_size = right_size = 0;
		if (q_buffer.size() < new_max_size) {
			q_buffer.resize(new_max_size);
			t_buffer.resize(new_max_size);
		}
	}
};

struct DPathData2 {
	int k, x1, y1, x2, y2, pre_k;
	explicit DPathData2() { }
	explicit DPathData2(const int j, const int ki, const int l, const int m, const int n, const int p) : k(j), x1(ki), y1(l), x2(m), y2(n), pre_k(p) { }
	void set(const int j, const int ki, const int l, const int m, const int n, const int p) {
		k = j;
		x1 = ki;
		y1 = l;
		x2 = m;
		y2 = n;
		pre_k = p;
	}
};

struct PathPoint {
	int x, y;
	explicit PathPoint() { }
	explicit PathPoint(const int i, const int j) : x(i), y(j) { }
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
	std::vector<char> query, target;
	std::vector<int> DynQ, DynT;
	std::vector<size_t> d_path_index;
	// maybe make d_path a deque rather than a vector, for growth
	// speed and (likely) reduced memory waste?
	std::vector<DPathData2> d_path;
	std::vector<PathPoint> aln_path;
    public:
	explicit DiffRunningData(const SW_Parameters& swp) : segment_size(swp.segment_size), DynQ(swp.row_size), DynT(swp.column_size) { }
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
