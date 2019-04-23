#include "dw.h"
#include <string.h>	// memcpy(), memset()
#include <vector>	// vector<>

#define GAP_ALN 4

static int CompareDPathData2(const void* const a, const void* const b) {
	const DPathData2* const d1((const DPathData2*)a);
	const DPathData2* const d2((const DPathData2*)b);
	return (d1->d != d2->d) ? (d1->d - d2->d) : (d1->k - d2->k);
}

static void fill_align(const std::string& query, const int q_offset, const std::string& target, const int t_offset, Alignment& align, std::vector<DPathData2>& d_path, std::vector<PathPoint>& aln_path, const int extend_forward, const size_t aln_path_max) {
	align.clear();
	align.aln_q_e = d_path.back().x2;
	align.aln_t_e = d_path.back().y2;
	align.dist = d_path.back().d;
	// get align path
	aln_path.clear();
	// no need to search for the first one
	const DPathData2* d_path_aux(&d_path[0] + d_path.size() - 1);
	aln_path.push_back(PathPoint(d_path_aux->x2, d_path_aux->y2));
	aln_path.push_back(PathPoint(d_path_aux->x1, d_path_aux->y1));
	DPathData2 seek(align.dist - 1, d_path.back().pre_k);
	for (; seek.d >= 0 && aln_path.size() < aln_path_max; --seek.d) {
		// there may be a better approach here than bsearch()
		d_path_aux = (const DPathData2*)bsearch(&seek, &d_path[0], d_path_aux - &d_path[0], sizeof(DPathData2), CompareDPathData2);
		aln_path.push_back(PathPoint(d_path_aux->x2, d_path_aux->y2));
		aln_path.push_back(PathPoint(d_path_aux->x1, d_path_aux->y1));
		seek.k = d_path_aux->pre_k;
	}
	std::vector<PathPoint>::const_reverse_iterator a(aln_path.rbegin());
	const std::vector<PathPoint>::const_reverse_iterator end_a(aln_path.rend());
	align.aln_q_s = a->x;
	align.aln_t_s = a->y;
	int current_x(a->x), current_y(a->y);
	for (++a; a != end_a; ++a) {
		const int new_x(a->x);
		const int new_y(a->y);
		if (current_x != new_x && current_y != new_y) {
			if (extend_forward) {
				align.q_aln_str.append(query, q_offset + current_x, new_x - current_x);
				align.t_aln_str.append(target, t_offset + current_y, new_y - current_y);
				current_x = new_x;
				current_y = new_y;
			} else {
				// don't know a clever way to append a reversed string
				for (; current_x < new_x; ++current_x) {
					align.q_aln_str += query[q_offset - current_x];
				}
				for (; current_y < new_y; ++current_y) {
					align.t_aln_str += target[t_offset - current_y];
				}
			}
		} else if (current_x != new_x) {
			align.t_aln_str.append(new_x - current_x, GAP_ALN);
			if (extend_forward) {
				align.q_aln_str.append(query, q_offset + current_x, new_x - current_x);
				current_x = new_x;
			} else {
				for (; current_x < new_x; ++current_x) {
					align.q_aln_str += query[q_offset - current_x];
				}
			}
		} else if (current_y != new_y) {
			align.q_aln_str.append(new_y - current_y, GAP_ALN);
			if (extend_forward) {
				align.t_aln_str.append(target, t_offset + current_y, new_y - current_y);
				current_y = new_y;
			} else {
				for (; current_y < new_y; ++current_y) {
					align.t_aln_str += target[t_offset - current_y];
				}
			}
		}
	}
}

static int Align(const std::string& query, const int q_offset, const int q_len, const std::string& target, const int t_offset, const int t_len, const int band_tolerance, Alignment& align, std::vector<int>& V, std::vector<int>& U, std::vector<DPathData2>& d_path, std::vector<PathPoint>& aln_path, const int extend_forward, const double error_rate) {
	const int k_offset(2 * error_rate * (q_len + t_len));
	const int band_size(band_tolerance * 2);
	d_path.clear();
	int best_m(-1), min_k(0), max_k(0);
	for (int d(0); d < k_offset && max_k - min_k <= band_size; ++d) {
		for (int k(min_k); k <= max_k; k += 2) {
			int x, pre_k;
			if (k == min_k || (k != max_k && V[k_offset + k - 1] < V[k_offset + k + 1])) {
				pre_k = k + 1;
				x = V[k_offset + k + 1];
			} else {
				pre_k = k - 1;
				x = V[k_offset + k - 1] + 1;
			}
			int y(x - k);
			const int x1(x), y1(y);
			if (extend_forward) {
				for (; x < q_len && y < t_len && query[q_offset + x] == target[t_offset + y]; ++x, ++y) { }
			} else {
				for (; x < q_len && y < t_len && query[q_offset - x] == target[t_offset - y]; ++x, ++y) { }
			}
			d_path.push_back(DPathData2(d, k, x1, y1, x, y, pre_k));
			if (x >= q_len || y >= t_len) {
				fill_align(query, q_offset, target, t_offset, align, d_path, aln_path, extend_forward, q_len + t_len + 1);
				return 1;
			}
			V[k_offset + k] = x;
			U[k_offset + k] = x + y;
			best_m = std::max(best_m, x + y);
		}
		// for banding
		int new_min_k(max_k);
		int new_max_k(min_k);
		for (int k(min_k); k <= max_k; k += 2) {
			if (U[k_offset + k] >= best_m - band_tolerance) {
				new_min_k = std::min(new_min_k, k);
				new_max_k = std::max(new_max_k, k);
			}
		}
		max_k = new_max_k + 1;
		min_k = new_min_k - 1;
	}
	return 0;
}

static void dw_in_one_direction(const std::string& query, const int q_offset, const std::string& target, const int t_offset, std::vector<int>& U, std::vector<int>& V, Alignment& align, std::vector<DPathData2>& d_path, std::vector<PathPoint>& aln_path, const SW_Parameters& swp, OutputStore& result, const int extend_forward, const double error_rate) {
	const int seg_size(swp.segment_size);
	int extend1(0), extend2(0);
	for (int not_at_end(1); not_at_end;) {
		// size left to extend
		int extend_size;
		if (extend_forward) {
			extend_size = std::min(query.size() - q_offset - extend1, target.size() - t_offset - extend2);
		} else {
			extend_size = std::min(q_offset - extend1, t_offset - extend2);
		}
		U.assign(swp.row_size, 0);
		V.assign(swp.column_size, 0);
		if (extend_size > seg_size + 100) {
			if (!Align(query, q_offset + (extend_forward ? extend1 : -extend1), seg_size, target, t_offset + (extend_forward ? extend2 : -extend2), seg_size, seg_size * 0.3, align, U, V, d_path, aln_path, extend_forward, error_rate)) {
				break;
			}
		} else {
			if (!Align(query, q_offset + (extend_forward ? extend1 : -extend1), extend_size, target, t_offset + (extend_forward ? extend2 : -extend2), extend_size, extend_size * 0.3, align, U, V, d_path, aln_path, extend_forward, error_rate)) {
				break;
			}
			not_at_end = 0;
		}
		int k, i(0), j(0), num_matches(0);
		for (k = align.q_aln_str.size() - 1; -1 < k && num_matches < 4; --k) {
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
			if (align.aln_q_e == i) {
				break;
			}
			extend1 += align.aln_q_e - i;
			extend2 += align.aln_t_e - j;
			++k;
		} else if (align.aln_q_e == 0) {
			break;
		} else {
			k = align.q_aln_str.size();
		}
		if (extend_forward) {
			result.right_store1 += align.q_aln_str.substr(0, k);
			result.right_store2 += align.t_aln_str.substr(0, k);
		} else {
			result.left_store1 += align.q_aln_str.substr(0, k);
			result.left_store2 += align.t_aln_str.substr(0, k);
		}
	}
}

static int dw(const std::string& query, const int query_start, const std::string& target, const int target_start, std::vector<int>& U, std::vector<int>& V, Alignment& align, std::vector<DPathData2>& d_path, std::vector<PathPoint>& aln_path, OutputStore& result, const SW_Parameters& swp, const double error_rate, const size_t min_aln_size) {
	result.clear();
	// reverse extend
	dw_in_one_direction(query, query_start - 1, target, target_start - 1, U, V, align, d_path, aln_path, swp, result, 0, error_rate);
	// forward extend
	dw_in_one_direction(query, query_start, target, target_start, U, V, align, d_path, aln_path, swp, result, 1, error_rate);
	// merge the results
	const char* const encode2char("ACGT-");
	int i, j, k;
	for (k = result.left_store1.size() - 1, i = 0, j = 0; -1 < k; --k) {
		int ch(result.left_store1[k]);
		if (ch < 0 || 4 < ch) {
			ERROR("Left1: Out of range 0-4: %d (%d, %s)", ch, k, result.left_store1.c_str());
		} else if (ch != 4) {	// not '-'
			++i;
		}
		result.out_store1 += encode2char[ch];
		ch = result.left_store2[k];
		if (ch < 0 || 4 < ch) {
			ERROR("Left2: Out of range 0-4: %d (%d, %s)", ch, k, result.left_store2.c_str());
		} else if (ch != 4) {	// not '-'
			++j;
		}
		result.out_store2 += encode2char[ch];
	}
	if (query_start < i) {
		ERROR("query_start %d, i %d, left_store_size1 %lu", query_start, i, result.left_store1.size());
	} else if (target_start < j) {
		ERROR("target_start %d, j %d, left_store_size2 %lu", target_start, j, result.left_store2.size());
	}
	result.query_start = query_start - i;
	result.target_start = target_start - j;
	for (k = 0, i = 0, j = 0; k < static_cast<int>(result.right_store1.size()); ++k) {
		int ch(result.right_store1[k]);
		if (ch < 0 || 4 < ch) {
			ERROR("Right1: Out of range 0-4: %d (%d, %s)", ch, k, result.right_store1.c_str());
		} else if (ch != 4) {	// not '-'
			++i;
		}
		result.out_store1 += encode2char[ch];
		ch = result.right_store2[k];
		if (ch < 0 || 4 < ch) {
			ERROR("Right2: Out of range 0-4: %d (%d, %s)", ch, k, result.right_store2.c_str());
		} else if (ch != 4) {	// not '-'
			++j;
		}
		result.out_store2 += encode2char[ch];
	}
	if (result.out_store1.size() < min_aln_size) {
		return 0;
	}
	result.query_end = query_start + i;
	result.target_end = target_start + j;
	for (size_t m(0); m < result.out_store1.size(); ++m) {
		if (result.out_store1[m] == result.out_store2[m]) {
			++result.mat;
			result.out_match_pattern += '|';
		} else if (result.out_store1[m] == '-') {
			++result.ins;
			result.out_match_pattern += '*';
		} else if (result.out_store2[m] == '-') {
			++result.del;
			result.out_match_pattern += '*';
		} else {
			++result.mis;
			result.out_match_pattern += '*';
		}
	}
	return 1;
}

bool GetAlignment(const std::string& query, const int query_start, const std::string& target, const int target_start, DiffRunningData& drd, M5Record& m5, const double error_rate, const int min_aln_size) {
	if (!dw(query, query_start, target, target_start, drd.DynQ, drd.DynT, drd.align, drd.d_path, drd.aln_path, drd.result, drd.swp, error_rate, min_aln_size)) {
		return 0;
	}
	const int consecutive_match_region_size(4);
	// trim starting end of alignment
	int qrb(0);	// q starting pads
	int trb(0);	// t starting pads
	int eit(0);	// matching run length
	size_t k;
	for (k = 0; k < drd.result.out_store1.size(); ++k) {
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
	const size_t start_aln_id(k - consecutive_match_region_size);
	// trim trailing end of alignment
	int qre(0);	// q ending pads
	int tre(0);	// t ending pads
	for (k = drd.result.out_store1.size() - 1, eit = 0; start_aln_id < k; --k) {
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
	const size_t aln_size(k + consecutive_match_region_size + 1 - start_aln_id);
	m5.m5qsize() = query.size();
	m5.m5qoff() = drd.result.query_start + qrb;
	m5.m5qend() = drd.result.query_end - qre;
	m5.m5qdir() = FWD;
	m5.m5ssize() = target.size();
	m5.m5soff() = drd.result.target_start + trb;
	m5.m5send() = drd.result.target_end - tre;
	m5.m5sdir() = FWD;
	m5.m5qaln() = drd.result.out_store1.substr(start_aln_id, aln_size);
	m5.m5saln() = drd.result.out_store2.substr(start_aln_id, aln_size);
	m5.m5pat() = drd.result.out_match_pattern.substr(start_aln_id, aln_size);
	return 1;
}
