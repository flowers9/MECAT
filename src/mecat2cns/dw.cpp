#include "dw.h"
#include <algorithm>	// copy(), fill()
#include <vector>	// vector<>

#define GAP_ALN 4

static void fill_align(const std::string& query, const int q_offset, const std::string& target, const int t_offset, Alignment& align, const std::vector<DPathData>& d_path, const DPathData* d_path_aux, const std::vector<DPathIndex>& d_path_index, int d, std::vector<PathPoint>& aln_path, const int extend_forward) {
	align.aln_q_e = d_path_aux->x2;
	align.aln_t_e = d_path_aux->y2;
	// get align path
	int aln_idx(-1);
	for (;;) {
		aln_path[++aln_idx].set(d_path_aux->x2, d_path_aux->y2);
		aln_path[++aln_idx].set(d_path_aux->x1, d_path_aux->y1);
		if (--d == -1) {
			break;
		}
		d_path_aux = &d_path[d_path_index[d].d_offset + (d_path_aux->pre_k - d_path_index[d].min_k) / 2];
	}
	// walk backwards along align path to fill in sequence with gaps
	align.reset();
	int current_x(aln_path[aln_idx].x);
	int current_y(aln_path[aln_idx].y);
	if (extend_forward) {
		const char* const query_p(query.data() + q_offset);
		const char* const target_p(target.data() + t_offset);
		for (--aln_idx; aln_idx != -1; --aln_idx) {
			const int new_x(aln_path[aln_idx].x);
			const int new_y(aln_path[aln_idx].y);
			const int dx(new_x - current_x);
			const int dy(new_y - current_y);
			if (dx && dy) {		// apparently, dx always equals dy in this case
				std::copy(query_p + current_x, query_p + new_x, &align.q_aln_str[align.size]);
				std::copy(target_p + current_y, target_p + new_y, &align.t_aln_str[align.size]);
				align.size += dx;
			} else if (dx) {
				std::copy(query_p + current_x, query_p + new_x, &align.q_aln_str[align.size]);
				std::fill(&align.t_aln_str[align.size], &align.t_aln_str[align.size] + dx, GAP_ALN);
				align.size += dx;
			} else if (dy) {
				std::fill(&align.q_aln_str[align.size], &align.q_aln_str[align.size] + dy, GAP_ALN);
				std::copy(target_p + current_y, target_p + new_y, &align.t_aln_str[align.size]);
				align.size += dy;
			}
			current_x = new_x;
			current_y = new_y;
		}
	} else {
		const char* const query_p(query.data() + q_offset + 1);
		const char* const target_p(target.data() + t_offset + 1);
		for (--aln_idx; aln_idx != -1; --aln_idx) {
			const int new_x(aln_path[aln_idx].x);
			const int new_y(aln_path[aln_idx].y);
			const int dx(new_x - current_x);
			const int dy(new_y - current_y);
			if (dx && dy) {		// apparently, dx always equals dy in this case
				align.size += dx;
				const int offset(align.q_aln_str.size() - align.size);
				std::copy(query_p - new_x, query_p - current_x, &align.q_aln_str[offset]);
				std::copy(target_p - new_y, target_p - current_y, &align.t_aln_str[offset]);
			} else if (dx) {
				align.size += dx;
				const int offset(align.q_aln_str.size() - align.size);
				std::copy(query_p - new_x, query_p - current_x, &align.q_aln_str[offset]);
				std::fill(&align.t_aln_str[offset], &align.t_aln_str[offset] + dx, GAP_ALN);
			} else if (dy) {
				align.size += dy;
				const int offset(align.q_aln_str.size() - align.size);
				std::fill(&align.q_aln_str[offset], &align.q_aln_str[offset] + dy, GAP_ALN);
				std::copy(target_p - new_y, target_p - current_y, &align.t_aln_str[offset]);
			}
			current_x = new_x;
			current_y = new_y;
		}
	}
}

static int Align(const int extend_size, const std::string& query, const int q_offset, const std::string& target, const int t_offset, Alignment& align, std::vector<int>& V, std::vector<int>& combined_match_length, std::vector<DPathData>& d_path, std::vector<DPathIndex>& d_path_index, std::vector<PathPoint>& aln_path, const int extend_forward, const double error_rate) {
	const int k_offset(extend_size * 4 * error_rate);
	const int band_tolerance(extend_size / 10 * 3);
	const int max_band_size(band_tolerance * 2 + 1);
	int d_path_idx(0), best_combined_match_length(0), min_k(0), max_k(0);
	V[k_offset + 1] = 0;	// initialize starting point
	for (int d(0); d != k_offset && max_k - min_k < max_band_size; ++d) {
		// starting point of each "d" set of entries
		d_path_index[d].set(d_path_idx, min_k);
		for (int k(min_k); k <= max_k; k += 2) {
			int q_pos, pre_k;
			if (k == min_k || (k != max_k && V[k_offset + k - 1] < V[k_offset + k + 1])) {
				pre_k = k + 1;
				q_pos = V[k_offset + k + 1];
			} else {
				pre_k = k - 1;
				q_pos = V[k_offset + k - 1] + 1;
			}
			int t_pos(q_pos - k);
			// start of exact match
			const int q_start(q_pos), t_start(t_pos);
			// find the other end of exact match
			if (extend_forward) {
				for (; q_pos < extend_size && t_pos < extend_size && query[q_offset + q_pos] == target[t_offset + t_pos]; ++q_pos, ++t_pos) { }
			} else {
				for (; q_pos < extend_size && t_pos < extend_size && query[q_offset - q_pos] == target[t_offset - t_pos]; ++q_pos, ++t_pos) { }
			}
			d_path[d_path_idx].set(q_start, t_start, q_pos, t_pos, pre_k);
			// see if we got as much as we can
			if (q_pos == extend_size || t_pos == extend_size) {
				fill_align(query, q_offset, target, t_offset, align, d_path, &d_path[d_path_idx], d_path_index, d, aln_path, extend_forward);
				return 1;
			}
			++d_path_idx;
			V[k_offset + k] = q_pos;
			combined_match_length[k_offset + k] = q_pos + t_pos;
			if (best_combined_match_length < q_pos + t_pos) {
				best_combined_match_length = q_pos + t_pos;
			}
		}
		// shift ends to one outside "good" band
		const int cutoff(best_combined_match_length - band_tolerance);
		for (; combined_match_length[k_offset + min_k] < cutoff; min_k += 2) { }
		--min_k;
		for (; combined_match_length[k_offset + max_k] < cutoff; max_k -= 2) { }
		++max_k;
	}
	return 0;
}

static void dw_in_one_direction(const std::string& query, const int q_offset, const std::string& target, const int t_offset, std::vector<int>& U, std::vector<int>& V, Alignment& align, std::vector<DPathData>& d_path, std::vector<DPathIndex>& d_path_index, std::vector<PathPoint>& aln_path, const int segment_size, OutputStore& result, const int extend_forward, const double error_rate) {
	const int q_extend_max(extend_forward ? query.size() - q_offset : q_offset + 1);
	const int t_extend_max(extend_forward ? target.size() - t_offset : t_offset + 1);
	// amount we've already extended
	int q_extend(0), t_extend(0), not_at_end(1);
	do {
		// size left to extend
		const int extend_size(std::min(q_extend_max - q_extend, t_extend_max - t_extend));
		if (extend_size <= segment_size + SEGMENT_BORDER) { // close enough to the end
			not_at_end = 0;
		}
		if (!Align(not_at_end ? segment_size : extend_size, query, q_offset + (extend_forward ? q_extend : -q_extend), target, t_offset + (extend_forward ? t_extend : -t_extend), align, U, V, d_path, d_path_index, aln_path, extend_forward, error_rate)) {
			return;
		}
		int k;
		if (not_at_end) {
			// go backwards until we find a good match (4 consecutive
			// matching basepairs), counting non-gap basepairs
			int q_bps(0), t_bps(0), num_matches(0);
			if (extend_forward) {
				for (k = align.size - 1; k != -1; --k) {
					const char qc(align.q_aln_str[k]);
					const char tc(align.t_aln_str[k]);
					if (qc != tc) {
						num_matches = 0;
						if (qc == GAP_ALN) {
							--q_bps;
						} else if (tc == GAP_ALN) {
							--t_bps;
						}
					} else if (++num_matches == 4) {
						break;
					}
				}
				q_bps += align.size - k;
				t_bps += align.size - k;
			} else {
				const int offset(align.q_aln_str.size() - align.size);
				for (k = offset; k != static_cast<int>(align.q_aln_str.size()); ++k) {
					const char qc(align.q_aln_str[k]);
					const char tc(align.t_aln_str[k]);
					if (qc != tc) {
						num_matches = 0;
						if (qc == GAP_ALN) {
							--q_bps;
						} else if (tc == GAP_ALN) {
							--t_bps;
						}
					} else if (++num_matches == 4) {
						++k;	// don't include final match
						break;
					}
				}
				q_bps += k - offset;
				t_bps += k - offset;
			}
			if (align.aln_q_e == q_bps) {	// no good match
				return;
			}
			// don't extend into the good match; keeping
			// good match for next alignment, maybe?
			q_extend += align.aln_q_e - q_bps;
			t_extend += align.aln_t_e - t_bps;
		} else if (align.aln_q_e == 0) {	// no good match
			return;
		} else if (extend_forward) {
			k = align.size;
		} else {
			k = align.q_aln_str.size() - align.size;
		}
		if (extend_forward) {
			const int offset(result.buffer_start + result.right_size);
			std::copy(&align.q_aln_str[0], &align.q_aln_str[0] + k, &result.q_buffer[offset]);
			std::copy(&align.t_aln_str[0], &align.t_aln_str[0] + k, &result.t_buffer[offset]);
			result.right_size += k;
		} else {
			result.left_size += align.q_aln_str.size() - k;
			const int offset(result.buffer_start - result.left_size);
			std::copy(&align.q_aln_str[k], &align.q_aln_str[0] + align.q_aln_str.size(), &result.q_buffer[offset]);
			std::copy(&align.t_aln_str[k], &align.t_aln_str[0] + align.q_aln_str.size(), &result.t_buffer[offset]);
		}
	} while (not_at_end);
}

static int gap_count(const std::vector<char>& buffer, int i, const int end_i) {
	int j(0);
	for (; i != end_i; ++i) {
		if (buffer[i] == GAP_ALN) {
			++j;
		}
	}
	return j;
}

static int dw(const std::string& query, const int query_start, const std::string& target, const int target_start, std::vector<int>& U, std::vector<int>& V, Alignment& align, std::vector<DPathData>& d_path, std::vector<DPathIndex>& d_path_index, std::vector<PathPoint>& aln_path, OutputStore& result, const int segment_size, const double error_rate, const int min_aln_size) {
	result.reset_buffer(query_start + target_start);
	// reverse extend (left side)
	dw_in_one_direction(query, query_start - 1, target, target_start - 1, U, V, align, d_path, d_path_index, aln_path, segment_size, result, 0, error_rate);
	// forward extend (right side)
	dw_in_one_direction(query, query_start, target, target_start, U, V, align, d_path, d_path_index, aln_path, segment_size, result, 1, error_rate);
	if (result.left_size + result.right_size < min_aln_size) {
		return 0;
	}
	const int left_start(result.buffer_start - result.left_size);
	const int right_end(result.buffer_start + result.right_size);
	result.query_start = query_start - result.left_size + gap_count(result.q_buffer, left_start, result.buffer_start);
	result.query_end = query_start + result.right_size - gap_count(result.q_buffer, result.buffer_start, right_end);
	result.target_start = target_start - result.left_size + gap_count(result.t_buffer, left_start, result.buffer_start);
	result.target_end = target_start + result.right_size - gap_count(result.t_buffer, result.buffer_start, right_end);
	return 1;
}

static void decode_sequence(std::string& out_seq, const std::vector<char>& in_seq, const size_t offset, const size_t size) {
	out_seq.resize(size);
	for (size_t i(0); i != size; ++i) {
		out_seq[i] = "ACGT-"[static_cast<int>(in_seq[offset + i])];
	}
}

int GetAlignment(const std::string& query, const int query_start, const std::string& target, const int target_start, DiffRunningData& drd, M5Record& m5, const double error_rate, const int min_aln_size) {
	if (!dw(query, query_start, target, target_start, drd.DynQ, drd.DynT, drd.align, drd.d_path, drd.d_path_index, drd.aln_path, drd.result, drd.segment_size, error_rate, min_aln_size)) {
		return 0;
	}
	// create return m5 record
	const OutputStore& result(drd.result);
	const int consecutive_match_region_size(4);
	// trim starting end of alignment
	int qrb(0);	// q starting basepair offset to good sequence
	int trb(0);	// t starting basepair offset to good sequence
	int eit(0);	// matching run length
	const int start_k(result.buffer_start - result.left_size);
	const int end_k(result.buffer_start + result.right_size);
	int k(start_k);
	for (; k != end_k; ++k) {
		const char qc(result.q_buffer[k]);
		const char tc(result.t_buffer[k]);
		if (qc != tc) {
			eit = 0;
			// we don't count gaps
			if (qc == GAP_ALN) {
				--qrb;
			} else if (tc == GAP_ALN) {
				--trb;
			}
		} else if (++eit == consecutive_match_region_size) {
			break;
		}
	}
	if (eit != consecutive_match_region_size) {	// no good match
		return 0;
	}
	// +1 for the ++k we skipped with the break
	const int start_aln_id(k + 1 - consecutive_match_region_size);
	qrb += start_aln_id - start_k;
	trb += start_aln_id - start_k;
	// trim trailing end of alignment
	int qre(0);	// q ending basepair offset to good sequence
	int tre(0);	// t ending basepair offset to good sequence
	eit = 0;	// still matching run length
	for (k = end_k - 1;; --k) {
		const char qc(result.q_buffer[k]);
		const char tc(result.t_buffer[k]);
		if (qc != tc) {
			eit = 0;
			// we don't count gaps
			if (qc == GAP_ALN) {
				--qre;
			} else if (tc == GAP_ALN) {
				--tre;
			}
		} else if (++eit == consecutive_match_region_size) {
			break;
		}
	}
	// --k we skipped in the break is countered by +1 for end-inclusion
	k += consecutive_match_region_size;
	qre += end_k - k;
	tre += end_k - k;
	m5.qoff = result.query_start + qrb;
	m5.qend = result.query_end - qre;
	m5.soff = result.target_start + trb;
	m5.send = result.target_end - tre;
	// +1 to make this end-inclusive
	const size_t aln_size(k - start_aln_id);
	decode_sequence(m5.qaln, result.q_buffer, start_aln_id, aln_size);
	decode_sequence(m5.saln, result.t_buffer, start_aln_id, aln_size);
	return 1;
}
