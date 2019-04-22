#ifndef DW_H
#define DW_H

#include <algorithm>
#include <string>	// string
#include <vector>	// vector<>

#include "../common/alignment.h"
#include "../common/defs.h"

struct SW_Parameters {
	idx_t segment_size;
	idx_t row_size;
	idx_t column_size;
	idx_t segment_aln_size;
	idx_t max_aln_size;
	SW_Parameters(const idx_t i, const idx_t j, const idx_t k, const idx_t l, const idx_t m) : segment_size(i), row_size(j), column_size(k), segment_aln_size(l), max_aln_size(m) { }
};

inline SW_Parameters get_sw_parameters_small() {
	// 1000 instead of 500 for "large"
	return SW_Parameters(500, 4096, 4096, 4096, 2 * MAX_SEQ_SIZE);
}

struct Alignment {
	int aln_str_size;
	int dist;
	int aln_q_s;
	int aln_q_e;
	int aln_t_s;
	int aln_t_e;
	char* q_aln_str;
	char* t_aln_str;
	void init() {
		aln_str_size = 0;
		aln_q_s = aln_q_e = 0;
		aln_t_s = aln_t_e = 0;
	}
	Alignment(const idx_t max_aln_size) {
		safe_malloc(q_aln_str, char, max_aln_size);
		safe_malloc(t_aln_str, char, max_aln_size);
	}
	~Alignment() {
		safe_free(q_aln_str);
		safe_free(t_aln_str);
	}
};

struct OutputStore {
	char* left_store1;
	char* left_store2;
	char* right_store1;
	char* right_store2;
	char* out_store1;
	char* out_store2;
	char* out_match_pattern;
	int left_store_size;
	int right_store_size;
	int out_store_size;
	int query_start, query_end;
	int target_start, target_end;
	int mat, mis, ins, del;
	double ident;
	OutputStore(const idx_t max_aln_size) {
		safe_malloc(left_store1, char, max_aln_size);
		safe_malloc(left_store2, char, max_aln_size);
		safe_malloc(right_store1, char, max_aln_size);
		safe_malloc(right_store2, char, max_aln_size);
		safe_malloc(out_store1, char, max_aln_size);
		safe_malloc(out_store2, char, max_aln_size);
		safe_malloc(out_match_pattern, char, max_aln_size);
	}
	~OutputStore() {
		safe_free(left_store1);
		safe_free(left_store2);
		safe_free(right_store1);
		safe_free(right_store2);
		safe_free(out_store1);
		safe_free(out_store2);
		safe_free(out_match_pattern);
	}
	void init() {
		left_store_size = right_store_size = out_store_size = 0;
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
	SW_Parameters swp;
	Alignment align;
	OutputStore result;
	std::string query, target;
	std::vector<int> DynQ, DynT;
	std::vector<DPathData2> d_path;
	std::vector<PathPoint> aln_path;
    public:
	explicit DiffRunningData(const SW_Parameters& swp_in) : swp(swp_in), align(swp_in.segment_aln_size), result(swp_in.max_aln_size) { }
	~DiffRunningData() { }
};

struct CandidateStartPosition {
	idx_t qoff;
	idx_t toff;
	idx_t tstart;
	idx_t tsize;
	idx_t tid;
	int left_q, left_t;
	int right_q, right_t;
	int num1, num2;
	int score;
	idx_t toff_in_aln;
	char chain;
};

bool GetAlignment(const char* query, int query_start, int query_size, const char* target, int target_start, int target_size, DiffRunningData& drd, M5Record& m5, double error_rate, int min_aln_size);

#endif  // DW_H
