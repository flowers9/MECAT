#ifndef DW_H
#define DW_H

#include <algorithm>
#include <string>	// string
#include <vector>	// vector<>

#include "../common/alignment.h"
#include "../common/defs.h"

struct SW_Parameters {
	int segment_size;
	int row_size;
	int column_size;
	SW_Parameters(const int i, const int j, const int k) : segment_size(i), row_size(j), column_size(k) { }
};

inline SW_Parameters get_sw_parameters_small() {
	// 1000 instead of 500 for "large"
	return SW_Parameters(500, 4096, 4096);
}

class Alignment {
    public:
	int dist;
	int aln_q_s, aln_q_e;
	int aln_t_s, aln_t_e;
	std::string q_aln_str;
	std::string t_aln_str;
    public:
	explicit Alignment() { }
	~Alignment() { }
	void clear() {
		q_aln_str.clear();
		t_aln_str.clear();
	}
};

class OutputStore {
    public:
	int query_start, query_end;
	int target_start, target_end;
	int mat, mis, ins, del;
	std::string left_store1, left_store2;
	std::string right_store1, right_store2;
	std::string out_store1, out_store2;
	std::string out_match_pattern;
    public:
	explicit OutputStore() { }
	~OutputStore() { }
	void clear() {
		left_store1.clear();
		left_store2.clear();
		right_store1.clear();
		right_store2.clear();
		out_store1.clear();
		out_store2.clear();
		out_match_pattern.clear();
	}
};

struct DPathData {
	int x1, y1, x2, y2, pre_k;
	explicit DPathData() { }
	explicit DPathData(const int i, const int j, const int k, const int l, const int m) : x1(i), y1(j), x2(k), y2(l), pre_k(m) { }
};

struct DPathData2 : public DPathData {
	int d, k;
	explicit DPathData2(const int i, const int j) : d(i), k(j) { }
	explicit DPathData2(const int i, const int j, const int k, const int l, const int m, const int n, const int p) : DPathData(k, l, m, n, p), d(i), k(j) { }
};

struct PathPoint {
	int x, y;
	explicit PathPoint(const int i, const int j) : x(i), y(j) { }
};

class DiffRunningData {
    public:
	const SW_Parameters swp;
	Alignment align;
	OutputStore result;
	std::string query, target;
	std::vector<int> DynQ, DynT;
	std::vector<DPathData2> d_path;
	std::vector<PathPoint> aln_path;
    public:
	explicit DiffRunningData(const SW_Parameters& swp_in) : swp(swp_in) { }
	~DiffRunningData() { }
};

bool GetAlignment(const std::string& query, int query_start, const std::string& target, int target_start, DiffRunningData& drd, M5Record& m5, double error_rate, int min_aln_size);

#endif  // DW_H
