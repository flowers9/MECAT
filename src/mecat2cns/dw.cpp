#include "dw.h"
#include <algorithm>	// copy(), fill()
#include <vector>	// vector<>

#define GAP_ALN 4

static int CompareDPathData2(const void* const a, const void* const b) {
	const DPathData2* const d1((const DPathData2*)a);
	const DPathData2* const d2((const DPathData2*)b);
	return (d1->d != d2->d) ? (d1->d - d2->d) : (d1->k - d2->k);
}

static void fill_align(const std::string& query, const int q_offset, const std::string& target, const int t_offset, Alignment& align, std::vector<DPathData2>& d_path, const size_t d_path_idx, std::vector<PathPoint>& aln_path, const int extend_forward, const size_t aln_path_max) {
	const DPathData2* d_path_aux(&d_path[d_path_idx - 1]);
	align.dist = d_path_aux->d;
	align.aln_q_e = d_path_aux->x2;
	align.aln_t_e = d_path_aux->y2;
	// get align path
	if (aln_path.size() < aln_path_max + 2) {
		aln_path.resize(aln_path_max + 2);
	}
	size_t aln_idx(-1);
	// no need to search for the first one
	aln_path[++aln_idx].set(d_path_aux->x2, d_path_aux->y2);
	aln_path[++aln_idx].set(d_path_aux->x1, d_path_aux->y1);
	DPathData2 seek(d_path_aux->d - 1, d_path_aux->pre_k);
	for (; -1 < seek.d && aln_idx < aln_path_max; --seek.d) {
		// there may be a better approach here than bsearch()
		d_path_aux = (const DPathData2*)bsearch(&seek, &d_path[0], d_path_aux - &d_path[0], sizeof(DPathData2), CompareDPathData2);
		aln_path[++aln_idx].set(d_path_aux->x2, d_path_aux->y2);
		aln_path[++aln_idx].set(d_path_aux->x1, d_path_aux->y1);
		seek.k = d_path_aux->pre_k;
	}
	int current_x(aln_path[aln_idx].x);
	int current_y(aln_path[aln_idx].y);
	align.reset(aln_path[0].x + aln_path[0].y);
	// +1 so [-dx, 0) becomes (-dx, 0], to match [0, dx)
	const char* const query_p(query.data() + q_offset + (extend_forward ? 0 : 1));
	const char* const target_p(target.data() + t_offset + (extend_forward ? 0 : 1));
	for (--aln_idx; aln_idx != size_t(-1); --aln_idx) {
		const int new_x(aln_path[aln_idx].x);
		const int new_y(aln_path[aln_idx].y);
		const int dx(new_x - current_x);
		const int dy(new_y - current_y);
		if (dx && dy) {
			// apparently, dx always equals dy in this case
			if (extend_forward) {
				std::copy(query_p + current_x, query_p + new_x, align.q_aln_str.begin() + align.size);
				std::copy(target_p + current_y, target_p + new_y, align.t_aln_str.begin() + align.size);
			} else {
				std::reverse_copy(query_p - new_x, query_p - current_x, align.q_aln_str.begin() + align.size);
				std::reverse_copy(target_p - new_y, target_p - current_y, align.t_aln_str.begin() + align.size);
			}
			align.size += dx;
		} else if (dx) {
			std::fill(align.t_aln_str.begin() + align.size, align.t_aln_str.begin() + align.size + dx, GAP_ALN);
			if (extend_forward) {
				std::copy(query_p + current_x, query_p + new_x, align.q_aln_str.begin() + align.size);
			} else {
				std::reverse_copy(query_p - new_x, query_p - current_x, align.q_aln_str.begin() + align.size);
			}
			align.size += dx;
		} else if (dy) {
			std::fill(align.q_aln_str.begin() + align.size, align.q_aln_str.begin() + align.size + dy, GAP_ALN);
			if (extend_forward) {
				std::copy(target_p + current_y, target_p + new_y, align.t_aln_str.begin() + align.size);
			} else {
				std::reverse_copy(target_p - new_y, target_p - current_y, align.t_aln_str.begin() + align.size);
			}
			align.size += dy;
		}
		current_x = new_x;
		current_y = new_y;
	}
}

static int Align(const int segment_size, const std::string& query, const int q_offset, const std::string& target, const int t_offset, Alignment& align, std::vector<int>& V, std::vector<int>& U, std::vector<DPathData2>& d_path, std::vector<PathPoint>& aln_path, const int extend_forward, const double error_rate) {
	const int k_offset(segment_size * 4 * error_rate);
	const int band_tolerance(segment_size / 10 * 3 + 1);
	const int max_band_size(band_tolerance * 2 - 1);
	size_t d_path_idx(0);
	int best_m(-1), min_k(0), max_k(0);
	for (int d(0); d < k_offset && max_k - min_k < max_band_size; ++d) {
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
				for (; x < segment_size && y < segment_size && query[q_offset + x] == target[t_offset + y]; ++x, ++y) { }
			} else {
				for (; x < segment_size && y < segment_size && query[q_offset - x] == target[t_offset - y]; ++x, ++y) { }
			}
			if (d_path_idx != d_path.size()) {
				d_path[d_path_idx].set(d, k, x1, y1, x, y, pre_k);
			} else {
				d_path.push_back(DPathData2(d, k, x1, y1, x, y, pre_k));
			}
			++d_path_idx;
			if (x == segment_size || y == segment_size) {
				fill_align(query, q_offset, target, t_offset, align, d_path, d_path_idx, aln_path, extend_forward, segment_size * 2);
				return 1;
			}
			V[k_offset + k] = x;
			U[k_offset + k] = x + y;
			if (best_m < x + y) {
				best_m = x + y;
			}
		}
		// for banding
		int new_min_k(max_k);
		int new_max_k(min_k);
		const int min_u(best_m - band_tolerance);
		for (int k(min_k); k <= max_k; k += 2) {
			if (min_u < U[k_offset + k]) {
				if (new_min_k > k) {
					new_min_k = k;
				}
				if (new_max_k < k) {
					new_max_k = k;
				}
			}
		}
		min_k = new_min_k - 1;
		max_k = new_max_k + 1;
	}
	return 0;
}

static void dw_in_one_direction(const std::string& query, const int q_offset, const std::string& target, const int t_offset, std::vector<int>& U, std::vector<int>& V, Alignment& align, std::vector<DPathData2>& d_path, std::vector<PathPoint>& aln_path, const int segment_size, OutputStore& result, const int extend_forward, const double error_rate) {
	int extend1(0), extend2(0);
	for (int not_at_end(1); not_at_end;) {
		// size left to extend
		int extend_size;
		if (extend_forward) {
			extend_size = std::min(query.size() - q_offset - extend1, target.size() - t_offset - extend2);
		} else {
			extend_size = std::min(q_offset - extend1, t_offset - extend2);
		}
		if (extend_size < segment_size + 101) {
			not_at_end = 0;
		}
		U.assign(U.size(), 0);
		V.assign(V.size(), 0);
		if (!Align(not_at_end ? segment_size : extend_size, query, q_offset + (extend_forward ? extend1 : -extend1), target, t_offset + (extend_forward ? extend2 : -extend2), align, U, V, d_path, aln_path, extend_forward, error_rate)) {
			break;
		}
		int k, i(0), j(0), num_matches(0);
		for (k = align.size - 1; -1 < k && num_matches < 4; --k) {
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
			k = align.size;
		}
		if (extend_forward) {
			std::copy(align.q_aln_str.begin(), align.q_aln_str.begin() + k, &result.q_buffer[result.buffer_start + result.right_size]);
			std::copy(align.t_aln_str.begin(), align.t_aln_str.begin() + k, &result.t_buffer[result.buffer_start + result.right_size]);
			result.right_size += k;
		} else {
			result.left_size += k;
			std::copy(align.q_aln_str.begin(), align.q_aln_str.begin() + k, &result.q_buffer[result.buffer_start - result.left_size]);
			std::copy(align.t_aln_str.begin(), align.t_aln_str.begin() + k, &result.t_buffer[result.buffer_start - result.left_size]);
		}

	}
}

static void count_basepairs(const OutputStore& result, int k, const int end_k, int i, int j) {
	for (; k != end_k; ++k) {
		int ch(result.q_buffer[k]);
		assert(-1 < ch && ch < 5);
		if (ch != GAP_ALN) {
			++i;
		}
		ch = result.t_buffer[k];
		assert(-1 < ch && ch < 5);
		if (ch != GAP_ALN) {
			++j;
		}
	}
}

static int dw(const std::string& query, const int query_start, const std::string& target, const int target_start, std::vector<int>& U, std::vector<int>& V, Alignment& align, std::vector<DPathData2>& d_path, std::vector<PathPoint>& aln_path, OutputStore& result, const int segment_size, const double error_rate, const int min_aln_size) {
	result.reset_buffer(query_start + target_start, query.size() + target.size());
	// reverse extend
	dw_in_one_direction(query, query_start - 1, target, target_start - 1, U, V, align, d_path, aln_path, segment_size, result, 0, error_rate);
	// forward extend
	dw_in_one_direction(query, query_start, target, target_start, U, V, align, d_path, aln_path, segment_size, result, 1, error_rate);
	if (result.left_size + result.right_size < min_aln_size) {
		return 0;
	}
	int i, j;
	// initialize i and j outside subroutine to avoid warning
	count_basepairs(result, result.buffer_start - result.left_size, result.buffer_start, i = 0, j = 0);
	result.query_start = query_start - i;
	result.target_start = target_start - j;
	count_basepairs(result, result.buffer_start, result.buffer_start + result.right_size, i = 0, j = 0);
	result.query_end = query_start + i;
	result.target_end = target_start + j;
	return 1;
}

static void decode_sequence(std::string& out_seq, const std::vector<char>& in_seq, const size_t offset, const size_t size) {
	out_seq.resize(size);
	for (size_t i(0); i != size; ++i) {
		out_seq[i] = "ACGT-"[static_cast<int>(in_seq[offset + i])];
	}
}

int GetAlignment(const std::string& query, const int query_start, const std::string& target, const int target_start, DiffRunningData& drd, M5Record& m5, const double error_rate, const int min_aln_size) {
	if (!dw(query, query_start, target, target_start, drd.DynQ, drd.DynT, drd.align, drd.d_path, drd.aln_path, drd.result, drd.segment_size, error_rate, min_aln_size)) {
		return 0;
	}
	const int consecutive_match_region_size(4);
	// trim starting end of alignment
	int qrb(0);	// q starting pads
	int trb(0);	// t starting pads
	int eit(0);	// matching run length
	int k;
	for (k = drd.result.buffer_start - drd.result.left_size; k != drd.result.buffer_start + drd.result.right_size; ++k) {
		const char qc(drd.result.q_buffer[k]);
		const char tc(drd.result.t_buffer[k]);
		if (qc != GAP_ALN) {
			++qrb;
		}
		if (tc != GAP_ALN) {
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
	eit = 0;	// still matching run length
	for (k = drd.result.buffer_start + drd.result.right_size - 1;; --k) {
		const char qc(drd.result.q_buffer[k]);
		const char tc(drd.result.t_buffer[k]);
		if (qc != GAP_ALN) {
			++qre;
		}
		if (tc != GAP_ALN) {
			++tre;
		}
		if (qc != tc) {
			eit = 0;
		} else if (++eit == consecutive_match_region_size) {
			break;
		}
	}
	qre -= consecutive_match_region_size;
	tre -= consecutive_match_region_size;
	m5.m5qsize() = query.size();
	m5.m5qoff() = drd.result.query_start + qrb;
	m5.m5qend() = drd.result.query_end - qre;
	m5.m5qdir() = FWD;
	m5.m5ssize() = target.size();
	m5.m5soff() = drd.result.target_start + trb;
	m5.m5send() = drd.result.target_end - tre;
	m5.m5sdir() = FWD;
	const size_t aln_size(k + consecutive_match_region_size - start_aln_id);
	decode_sequence(m5.m5qaln(), drd.result.q_buffer, start_aln_id, aln_size);
	decode_sequence(m5.m5saln(), drd.result.t_buffer, start_aln_id, aln_size);
	return 1;
}
