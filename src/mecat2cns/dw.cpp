#include "dw.h"
#include <string.h>	// memcpy(), memset()
#include <vector>	// vector<>

#define GAP_ALN 4

static int CompareDPathData2(const void* const a, const void* const b) {
	const DPathData2* const d1((const DPathData2*)a);
	const DPathData2* const d2((const DPathData2*)b);
	return (d1->d != d2->d) ? (d1->d - d2->d) : (d1->k - d2->k);
}

static void fill_align(const char* const query, const char* const target, Alignment& align, std::vector<DPathData2>& d_path, std::vector<PathPoint>& aln_path, const int right_extend, const size_t aln_path_max) {
	align.init();
	align.aln_q_e = d_path.back().x2;
	align.aln_t_e = d_path.back().y2;
	align.dist = d_path.back().d;
	align.aln_str_size = (align.aln_q_e + align.aln_t_e + align.dist) / 2;
	// get align path
	aln_path.clear();
	DPathData2 seek(align.dist, d_path.back().k);
	for (; seek.d >= 0 && aln_path.size() < aln_path_max; --seek.d) {
		// there may be a better approach here than bsearch()
		const DPathData2* const d_path_aux((const DPathData2*)bsearch(&seek, &d_path[0], d_path.size(), sizeof(DPathData2), CompareDPathData2));
		aln_path.push_back(PathPoint(d_path_aux->x2, d_path_aux->y2));
		aln_path.push_back(PathPoint(d_path_aux->x1, d_path_aux->y1));
		seek.k = d_path_aux->pre_k;
	}
	std::vector<PathPoint>::const_reverse_iterator a(aln_path.rbegin());
	const std::vector<PathPoint>::const_reverse_iterator end_a(aln_path.rend());
	int current_x(a->x);
	int current_y(a->y);
	align.aln_q_s = current_x;
	align.aln_t_s = current_y;
	int aln_pos(0);
	// starting increment is safe as we're guaranteed two entries at least
	for (++a; a != end_a; ++a) {
		const int new_x(a->x);
		const int new_y(a->y);
		const int dx(new_x - current_x);
		const int dy(new_y - current_y);
		if (dx == 0 && dy == 0) {
			continue;
		} else if (dx == 0 && dy != 0) {
			if (right_extend) {
				for (int i(0); i < dy; ++i) {
					align.q_aln_str[aln_pos + i] = GAP_ALN;
					align.t_aln_str[aln_pos + i] = target[current_y + i];
				}
			} else {
				for (int i(0); i < dy; ++i) {
					align.q_aln_str[aln_pos + i] = GAP_ALN;
					align.t_aln_str[aln_pos + i] = target[-(current_y + i)];
				}
			}
			aln_pos += dy;
		} else if (dx != 0 && dy == 0) {
			if (right_extend) {
				for (int i(0); i < dx; ++i) {
					align.q_aln_str[aln_pos + i] = query[current_x + i];
					align.t_aln_str[aln_pos + i] = GAP_ALN;
				}
			} else {
				for (int i(0); i < dx; ++i) {
					align.q_aln_str[aln_pos + i] = query[-(current_x + i)];
					align.t_aln_str[aln_pos + i] = GAP_ALN;
				}
			}
			aln_pos += dx;
		} else {
			if (right_extend) {
				for (int i(0); i < dx; ++i) {
					align.q_aln_str[aln_pos + i] = query[current_x + i];
				}
				for (int i(0); i < dy; ++i) {
					align.t_aln_str[aln_pos + i] = target[current_y + i];
				}
			} else {
				for (int i(0); i < dx; ++i) {
					align.q_aln_str[aln_pos + i] = query[-(current_x + i)];
				}
				for (int i(0); i < dy; ++i) {
					align.t_aln_str[aln_pos + i] = target[-(current_y + i)];
				}
			}
			aln_pos += dy;
		}
		current_x = new_x;
		current_y = new_y;
	}
	align.aln_str_size = aln_pos;
}

static int Align(const char* const query, const int q_len, const char* const target, const int t_len, const int band_tolerance, Alignment& align, std::vector<int>& V, std::vector<int>& U, std::vector<DPathData2>& d_path, std::vector<PathPoint>& aln_path, const int right_extend, const double error_rate) {
	const int k_offset(2 * error_rate * (q_len + t_len));
	const int band_size(band_tolerance * 2);
	d_path.clear();
	int best_m(-1), min_k(0), max_k(0);
	for (int d(0); d < k_offset && max_k - min_k <= band_size; ++d) {
		for (int k(min_k); k <= max_k; k += 2) {
			int x, pre_k;
			if (k == min_k || (k != max_k && V[k - 1 + k_offset] < V[k + 1 + k_offset])) {
				pre_k = k + 1;
				x = V[k + 1 + k_offset];
			} else {
				pre_k = k - 1;
				x = V[k - 1 + k_offset] + 1;
			}
			int y(x - k);
			const int x1(x), y1(y);
			if (right_extend) {
				while (x < q_len && y < t_len && query[x] == target[y]) {
					++x;
					++y;
				}
			} else {
				while (x < q_len && y < t_len && query[-x] == target[-y]) {
					++x;
					++y;
				}
			}
			d_path.push_back(DPathData2(d, k, x1, y1, x, y, pre_k));
			if (x >= q_len || y >= t_len) {
				fill_align(query, target, align, d_path, aln_path, right_extend, q_len + t_len + 1);
				return 1;
			}
			V[k + k_offset] = x;
			U[k + k_offset] = x + y;
			best_m = std::max(best_m, x + y);
		}
		// for banding
		int new_min_k(max_k);
		int new_max_k(min_k);
		for (int k2(min_k); k2 <= max_k; k2 += 2) {
			if (U[k2 + k_offset] >= best_m - band_tolerance) {
				new_min_k = std::min(new_min_k, k2);
				new_max_k = std::max(new_max_k, k2);
			}
		}
		max_k = new_max_k + 1;
		min_k = new_min_k - 1;
	}
	return 0;
}

static void dw_in_one_direction(const char* const query, const int query_size, const char* const target, const int target_size, std::vector<int>& U, std::vector<int>& V, Alignment& align, std::vector<DPathData2>& d_path, std::vector<PathPoint>& aln_path, const SW_Parameters& swp, OutputStore& result, const int right_extend, const double error_rate) {
	const idx_t ALN_SIZE(swp.segment_size);
	const idx_t U_SIZE(swp.row_size);
	const idx_t V_SIZE(swp.column_size);
	int extend_size(std::min(query_size, target_size));	// size left to extend
	int seg_size(ALN_SIZE);
	int extend1(0), extend2(0);
	for (int not_at_end(1); not_at_end;) {
		if (extend_size <= ALN_SIZE + 100) {
			seg_size = extend_size;
			not_at_end = 0;
		}
		const char* seq1;
		const char* seq2;
		if (right_extend) {
			seq1 = query + extend1;
			seq2 = target + extend2;
		} else {
			seq1 = query - extend1;
			seq2 = target - extend2;
		}
		U.assign(U_SIZE, 0);
		V.assign(V_SIZE, 0);
		if (!Align(seq1, seg_size, seq2, seg_size, seg_size * 0.3, align, U, V, d_path, aln_path, right_extend, error_rate)) {
			break;
		}
		int i(0), j(0), k, num_matches(0);
		for (k = align.aln_str_size - 1; k > -1 && num_matches < 4; --k) {
			if (align.q_aln_str[k] != GAP_ALN) {
				++i;
			}
			if (align.t_aln_str[k] != GAP_ALN) {
				++j;
			}
			if (align.q_aln_str[k] == align.t_aln_str[k]) {
				++num_matches;
			} else {
				num_matches = 0;
			}
		}
		if (not_at_end) {
			++k;
			i += ALN_SIZE - align.aln_q_e;
			if (i == ALN_SIZE) {
				break;
			}
			j += ALN_SIZE - align.aln_t_e;
			extend1 += ALN_SIZE - i;
			extend2 += ALN_SIZE - j;
		} else {
			if (align.aln_q_e == 0) {
				break;
			}
			extend1 += align.aln_q_e;
			extend2 += align.aln_t_e;
			k = align.aln_str_size;
		}
		if (right_extend) {
			memcpy(result.right_store1 + result.right_store_size, align.q_aln_str, k);
			memcpy(result.right_store2 + result.right_store_size, align.t_aln_str, k);
			result.right_store_size += k;
		} else {
			memcpy(result.left_store1 + result.left_store_size, align.q_aln_str, k);
			memcpy(result.left_store2 + result.left_store_size, align.t_aln_str, k);
			result.left_store_size += k;
		}
		extend_size = std::min(query_size - extend1, target_size - extend2);
	}
}

static int dw(const char* query, const int query_size, const int query_start, const char* target, const int target_size, const int target_start, std::vector<int>& U, std::vector<int>& V, Alignment& align, std::vector<DPathData2>& d_path, std::vector<PathPoint>& aln_path, OutputStore& result, const SW_Parameters& swp, const double error_rate, const int min_aln_size) {
	result.init();
	align.init();
	// left extend
	dw_in_one_direction(query + query_start - 1, query_start, target + target_start - 1, target_start, U, V, align, d_path, aln_path, swp, result, 0, error_rate);
	align.init();
	// right extend
	dw_in_one_direction(query + query_start, query_size - query_start, target + target_start, target_size - target_start, U, V, align, d_path, aln_path, swp, result, 1, error_rate);
	// merge the results
	int i, j, k, idx(0);
	const char* const encode2char("ACGT-");
	for (k = result.left_store_size - 1, i = 0, j = 0; -1 < k; --k, ++idx) {
		int ch(result.left_store1[k]);
		if (ch < 0 || 4 < ch) {
			ERROR("Left1: Out of range 0-4: %d (%d, %d, %d)", ch, k, idx, result.left_store_size);
		}
		ch = encode2char[ch];
		result.out_store1[idx] = ch;
		if (ch != '-') {
			++i;
		}
		ch = result.left_store2[k];
		if (ch < 0 || 4 < ch) {
			ERROR("Left2: Out of range 0-4: %d (%d, %d, %d)", ch, k, idx, result.left_store_size);
		}
		ch = encode2char[ch];
		result.out_store2[idx] = ch;
		if (ch != '-') {
			++j;
		}
	}
	if (query_start < i) {
		ERROR("query_start %d, i %d", query_start, i);
	} else if (target_start < j) {
		ERROR("target_start %d, j %d", target_start, j);
	}
	result.query_start = query_start - i;
	result.target_start = target_start - j;
	for (k = 0, i = 0, j = 0; k < result.right_store_size; ++k, ++idx) {
		int ch(result.right_store1[k]);
		if (ch < 0 || 4 < ch) {
			ERROR("Right1: Out of range 0-4: %d (%d, %d, %d)", ch, k, idx, result.right_store_size);
		}
		ch = encode2char[ch];
		result.out_store1[idx] = ch;
		if (ch != '-') {
			++i;
		}
		ch = result.right_store2[k];
		if (ch < 0 || 4 < ch) {
			ERROR("Right2: Out of range 0-4: %d (%d, %d, %d)", ch, k, idx, result.right_store_size);
		}
		ch = encode2char[ch];
		result.out_store2[idx] = ch;
		if (ch != '-') {
			++j;
		}
	}
	if (idx < min_aln_size) {
		return 0;
	}
	result.out_store_size = idx;
	result.query_end = query_start + i;
	result.target_end = target_start + j;
	int mat(0), mis(0), ins(0), del(0);
	for (j = 0; j < result.out_store_size; ++j) {
		if (result.out_store1[j] == result.out_store2[j]) {
			++mat;
			result.out_match_pattern[j] = '|';
		} else if (result.out_store1[j] == '-') {
			++ins;
			result.out_match_pattern[j] = '*';
		} else if (result.out_store2[j] == '-') {
			++del;
			result.out_match_pattern[j] = '*';
		} else {
			++mis;
			result.out_match_pattern[j] = '*';
		}
	}
	result.out_store1[result.out_store_size] = 0;
	result.out_store2[result.out_store_size] = 0;
	result.out_match_pattern[result.out_store_size] = 0;
	result.mat = mat;
	result.mis = mis;
	result.ins = ins;
	result.del = del;
	result.ident = double(100) * mat / result.out_store_size;
	return 1;
}

bool GetAlignment(const char* const query, const int query_start, const int query_size, const char* const target, const int target_start, const int target_size, DiffRunningData& drd, M5Record& m5, const double error_rate, const int min_aln_size) {
	if (!dw(query, query_size, query_start, target, target_size, target_start, drd.DynQ, drd.DynT, drd.align, drd.d_path, drd.aln_path, drd.result, drd.swp, error_rate, min_aln_size)) {
		return 0;
	}
	const int consecutive_match_region_size(4);
	// trim starting end of alignment
	int qrb(0);	// q starting pads
	int trb(0);	// t starting pads
	int eit(0);	// matching run length
	int k;
	for (k = 0; k < drd.result.out_store_size; ++k) {
		const char qc(drd.result.out_store1[k]);
		const char tc(drd.result.out_store2[k]);
		if (qc != '-') {
			++qrb;
		}
		if (tc != '-') {
			++trb;
		}
		if (qc != tc) {
			eit = 0;
		} else if (++eit == consecutive_match_region_size) {
			++k;
			break;
		}
	}
	if (eit < consecutive_match_region_size) {	// no good match
		return 0;
	}
	qrb -= consecutive_match_region_size;
	trb -= consecutive_match_region_size;
	const int start_aln_id(k - consecutive_match_region_size);
	// trim trailing end of alignment
	int qre(0);	// q ending pads
	int tre(0);	// t ending pads
	for (k = drd.result.out_store_size - 1, eit = 0; start_aln_id < k; --k) {
		const char qc(drd.result.out_store1[k]);
		const char tc(drd.result.out_store2[k]);
		if (qc != '-') {
			++qre;
		}
		if (tc != '-') {
			++tre;
		}
		if (qc != tc) {
			eit = 0;
		} else if (++eit == consecutive_match_region_size) {
			--k;
			break;
		}
	}
	qre -= consecutive_match_region_size;
	tre -= consecutive_match_region_size;
	const int end_aln_id(k + consecutive_match_region_size + 1);
	m5qsize(m5) = query_size;
	m5qoff(m5) = drd.result.query_start + qrb;
	m5qend(m5) = drd.result.query_end - qre;
	m5qdir(m5) = FWD;
	m5ssize(m5) = target_size;
	m5soff(m5) = drd.result.target_start + trb;
	m5send(m5) = drd.result.target_end - tre;
	m5sdir(m5) = FWD;
	const int aln_size(end_aln_id - start_aln_id);
	m5.m5qaln().assign(drd.result.out_store1 + start_aln_id, aln_size);
	m5.m5saln().assign(drd.result.out_store2 + start_aln_id, aln_size);
	m5.m5pat().assign(drd.result.out_match_pattern + start_aln_id, aln_size);
	return 1;
}
