#include "mecat_correction.h"

#include "MECAT_AlnGraphBoost.H"

using namespace ns_banded_sw;

namespace ns_meap_cns {

#define FMAT 1
#define FDEL 2
#define FINS 4
#define UNDS 8

// returns type of coverage present

inline uint1 identify_one_consensus_item(const CnsTableItem& cns_item, const int min_cov) {
	const int cov(ceil((cns_item.mat_cnt + cns_item.ins_cnt) * 0.8));
	uint1 ident;
	if (cns_item.mat_cnt >= cov) {		// coverage is 80% or more matches
		ident = FMAT;
	} else if (cns_item.ins_cnt >= cov) {	// coverage is 80% or more inserts
		ident = FINS;
	} else {				// neither of the above
		ident = UNDS;
	}
	if (2 * cns_item.del_cnt >= cov) {	// deletes are 40% or more than the coverage
		ident |= FDEL;
	}
	return ident;
}

struct CompareOverlapByOverlapSize
{
	bool operator()(const Overlap& a ,const Overlap& b)
	{
		const idx_t ovlp_a = std::max(a.qend - a.qoff, a.send - a.soff);
		const idx_t ovlp_b = std::max(b.qend - b.qoff, b.send - b.soff);
		return ovlp_a > ovlp_b;
	}
};

void meap_add_one_aln(const std::string& qaln, const std::string& saln, idx_t start_soff, CnsTableItem* const cns_table, const char* const org_seq) {
	r_assert(qaln.size() == saln.size());
	const idx_t aln_size(qaln.size());
	for (idx_t i(0); i < aln_size; ) {
		const char q(qaln[i]);
		const char s(saln[i]);
		if (q == '-' && s == '-') {	// skip
			++i;
		} else if (q == s) {		// match
			++cns_table[start_soff].mat_cnt;
			cns_table[start_soff].base = s;
			++start_soff;
			++i;
		} else if (q == '-') {		// insert
			++cns_table[start_soff].ins_cnt;
			++start_soff;
			++i;
		} else {			// delete
			r_assert(s == '-');
			for (++i; i < aln_size && saln[i] == '-'; ++i) { }
			++cns_table[start_soff - 1].del_cnt;
		}
	}
}

void meap_cns_one_indel(const int sb, const int se, CnsAlns& cns_vec, const int min_cov, std::string& aux_qstr, std::string& aux_tstr, std::string& cns) {
	AlnGraphBoost ag(se - sb + 1);
	int sb_out;
	for (CnsAln* a(cns_vec.begin()); a != cns_vec.end(); ++a) {
		if (a->retrieve_aln_subseqs(sb, se, aux_qstr, aux_tstr, sb_out)) {
			ag.addAln(aux_qstr, aux_tstr, sb_out - sb + 1);
		}
	}		
	ag.mergeNodes();
	ag.consensus(min_cov * 0.4, cns);
}

void meap_consensus_one_segment(const CnsTableItem* const cns_list, const int cns_list_size, uint1* const cns_id_vec, const int start_soff, CnsAlns& cns_vec, std::string& aux_qstr, std::string& aux_tstr, std::string& target, const int min_cov) {
	// get types of coverage
	for (int i(0); i < cns_list_size; ++i) {
		cns_id_vec[i] = identify_one_consensus_item(cns_list[i], min_cov);
	}
	std::string cns;
	target.clear();
	const uint1 unds_or_fdel(UNDS | FDEL);	// questionable coverage types
	int i(0); 
	// advance to matching coverage
	for (; i < cns_list_size && !(cns_id_vec[i] & FMAT); ++i) { }
	while (i < cns_list_size) {
		target.push_back(cns_list[i].base);
		const int start(i);
		// advance to next matching coverage
		for (++i; i < cns_list_size && !(cns_id_vec[i] & FMAT); ++i) { }
		int need_refinement(0);
		// check to see if anything in-between has questionable coverage
		for (int k(start); k < i; ++k) {
			if (cns_id_vec[k] & unds_or_fdel) {
				need_refinement = 1;
				break;
			}
		}
		if (need_refinement) {
			meap_cns_one_indel(start + start_soff, i + start_soff, cns_vec, cns_list[start].mat_cnt + cns_list[start].ins_cnt, aux_qstr, aux_tstr, cns);
			// trim first and last as they have good coverage
			if (cns.size() > 2) {
				target.append(cns.data() + 1, cns.size() - 2);
			}
		}
	}
}

struct CmpMappingRangeBySoff {
	bool operator()(const MappingRange& a, const MappingRange& b) {
		if (a.start != b.start) {
			return a.start < b.start;
		} else {
			return b.end < a.end;
		}
	}
};

void get_effective_ranges(std::vector<MappingRange>& mranges, std::vector<MappingRange>& eranges, const int read_size, const int min_size) {
	eranges.clear();
	if (mranges.size() == 0) {
		return;
	}
	// see if any ranges are effectively the whole read
	std::vector<MappingRange>::const_iterator a(mranges.begin());
	const std::vector<MappingRange>::const_iterator end_a(mranges.end());
	for (; a != end_a; ++a) {
		if (a->start <= 500 && read_size - a->end <= 500) {
			eranges.push_back(MappingRange(0, read_size));
			return;
		}
	}
	std::sort(mranges.begin(), mranges.end(), CmpMappingRangeBySoff());
	const int nr(mranges.size());
	// -1 so we can use > instead of >=
	const int min_size_95(ceil(min_size * 0.95) - 1);
	int i(0), left(mranges[0].start);
	for (;;) {
		const int right(mranges[i].end);
		// include ranges we completely overlap
		for (++i; i < nr && mranges[i].end <= right; ++i) { }
		if (i == nr) {
			if (right - left > min_size_95) {
				eranges.push_back(MappingRange(left, right));
			}
			break;
		}
		// if the next non-contained range overlaps by 1k or more,
		// treat as part of this range and keep extending
		if (right - mranges[i].start < 1000) {
			// truncate current range at overlap start
			const int eff_right(std::min(right, mranges[i].start));
			if (eff_right - left > min_size_95) {
				eranges.push_back(MappingRange(left, eff_right));
			}
			// put some space between the effective ranges
			left = std::max(right, mranges[i].start);
		}
	}
}

// breaks output sequence into chunks of no more than MaxSeqSize (if needed);
// split chunks will have an OvlpSize overlap

void
output_cns_result(std::vector<CnsResult>& cns_results,
				  CnsResult& cr,
				  const idx_t beg,
				  const idx_t end,
				  const std::string& cns_seq) 
{
	const size_t MaxSeqSize = 60000;
	const size_t OvlpSize = 10000;
	// BlkSize must be >= OvlpSize
	const size_t BlkSize = MaxSeqSize - OvlpSize - 1000;
	
	const size_t size = cns_seq.size();
	if (size <= MaxSeqSize) {
		cr.range[0] = beg;
		cr.range[1] = end;
		cr.seq = cns_seq;
		cns_results.push_back(cr);
	} else {
		const size_t cutoff = size - OvlpSize - 1000;
		size_t L = 0, R;
		do {
			R = L + BlkSize;
			if (R >= cutoff) {
				R = size;
			}
			cr.range[0] = L + beg;
			cr.range[1] = R < size && R + beg < static_cast<size_t>(end) ? R + beg : end;
			cr.seq = cns_seq.substr(L, R - L);
			cns_results.push_back(cr);
			L = R - OvlpSize;
		} while (R < size);
	}
}

inline bool
check_ovlp_mapping_range(const int qb, const int qe, const int qs,
						 const int sb, const int se, const int ss,
						 double ratio)
{
	const int oq = qe - qb;
	const int qqs = qs * ratio;
	const int os = se - sb;
	const int qss = ss * ratio;
	return oq >= qqs || os >= qss;
}

// look for areas of high coverage of about min_size or more,
// improve them and stick on the results pile

void consensus_worker(const CnsTableItem* const cns_table, uint1* const id_list, CnsAlns& cns_vec, std::string& aux_qstr, std::string& aux_tstr, const std::vector<MappingRange>& eranges, const int min_cov, const int min_size, const int read_id, std::vector<CnsResult>& cns_results) {
	CnsResult cns_result;
	cns_result.id = read_id;
	const idx_t min_size_95(ceil(0.95 * min_size));
	std::string cns_seq;
	std::vector<MappingRange>::const_iterator a(eranges.begin());
	const std::vector<MappingRange>::const_iterator end_a(eranges.end());
	for (; a != end_a; ++a) {
		const int end_i(a->end);
		for (idx_t i(a->start); i < end_i;) {
			// find start of next high coverage area
			for (; i < end_i && cns_table[i].mat_cnt + cns_table[i].ins_cnt < min_cov; ++i) { }
			const idx_t start(i);
			// find end of high coverage area
			for (++i; i < end_i && cns_table[i].mat_cnt + cns_table[i].ins_cnt >= min_cov; ++i) { }
			if (i - start >= min_size_95) {
				meap_consensus_one_segment(cns_table + start, i - start, id_list, start, cns_vec, aux_qstr, aux_tstr, cns_seq, min_cov);
				if (cns_seq.size() >= static_cast<size_t>(min_size)) {
					output_cns_result(cns_results, cns_result, start, i, cns_seq);
				}
			}
		}
	}
}

static void decode_and_append_sequence(std::string& s, const char* const seq, idx_t i, const idx_t end_i) {
	s.reserve(s.size() + end_i - i);
	for (; i < end_i; ++i) {
		s += "ACGT"[static_cast<int>(seq[i])];
	}
}

// same as consensus_worker, but produces entire read as one entry;
// uncorrected sections are just copied as is;

void consensus_worker_one_read(const CnsTableItem* const cns_table, uint1* const id_list, CnsAlns& cns_vec, std::string& aux_qstr, std::string& aux_tstr, const std::vector<MappingRange>& eranges, const int min_cov, const int min_size, const int read_id, const std::vector<char>& tstr, std::vector<CnsResult>& cns_results) {
	CnsResult cns_result;
	cns_result.id = read_id;
	const idx_t min_size_95(std::max<int>(ceil(0.95 * min_size), 1) - 1);
	std::string cns_seq;
	std::vector<MappingRange>::const_iterator a(eranges.begin());
	const std::vector<MappingRange>::const_iterator end_a(eranges.end());
	std::vector<MappingRange>::const_iterator last_a(eranges.end());
	for (; a != end_a; last_a = a++) {
		if (last_a != end_a) {		// add in-between range to cns_result
			decode_and_append_sequence(cns_result.seq, tstr.data(), last_a->end, a->start);
		} else if (a->start > 0) {	// add beginning of read
			decode_and_append_sequence(cns_result.seq, tstr.data(), 0, a->start);
		}
		const int begin_i(a->start - 1);
		const int end_i(a->end);
		for (idx_t i(begin_i); i < end_i;) {
			// find start of next high coverage area
			const idx_t last_end(i != begin_i ? i : a->start);
			for (++i; i < end_i && cns_table[i].mat_cnt + cns_table[i].ins_cnt < min_cov; ++i) { }
			// add low coverage area as-is
			decode_and_append_sequence(cns_result.seq, tstr.data(), last_end, i);
			if (i == end_i) {
				break;
			}
			const idx_t start(i);
			// find end of high coverage area
			for (++i; i < end_i && cns_table[i].mat_cnt + cns_table[i].ins_cnt >= min_cov; ++i) { }
			if (i - start > min_size_95) {
				meap_consensus_one_segment(cns_table + start, i - start, id_list, start, cns_vec, aux_qstr, aux_tstr, cns_seq, min_cov);
				if (cns_seq.size() >= static_cast<size_t>(min_size)) {
					// add corrected sequence
					cns_result.seq += cns_seq;
					continue;
				}
			}
			// add uncorrected sequence
			decode_and_append_sequence(cns_result.seq, tstr.data(), start, i);
		}
	}
	// add end of read
	if (last_a != end_a) {
		decode_and_append_sequence(cns_result.seq, tstr.data(), last_a->end, tstr.size());
	}
	cns_result.range[0] = 0;
	cns_result.range[1] = cns_result.seq.size();
	cns_results.push_back(cns_result);
}

void consensus_one_read_m4_pacbio(ConsensusThreadData& ctd, ConsensusPerThreadData &pctd, const idx_t read_id, const idx_t sid, const idx_t eid) {
	PackedDB& reads(ctd.reads);
	ExtensionCandidate* overlaps((ExtensionCandidate*)pctd.candidates);
	DiffRunningData* const drd_s(pctd.drd_s);
	DiffRunningData* drd(NULL);
	M5Record& m5(pctd.m5);
	CnsAlns& cns_vec(pctd.cns_alns);
	std::vector<CnsResult>& cns_results(pctd.cns_results);
	const idx_t read_size(overlaps[read_id].ssize);
	std::vector<char>& qstr(pctd.query);
	std::vector<char>& tstr(pctd.target);
	tstr.resize(read_size);
	reads.GetSequence(read_id, true, tstr.data(), read_size);
	std::string& nqstr(pctd.qaln);
	std::string& ntstr(pctd.saln);
	const int min_align_size(ctd.rco.min_align_size);
	const int max_added(60);
	CnsTableItem* cns_table(pctd.cns_table);
	std::for_each(cns_table, cns_table + read_size, CnsTableItemCleaner());
	cns_vec.clear();
	const idx_t L(sid);
	const idx_t R(eid - sid <= max_added ? eid : L + max_added);
	if (eid - sid > max_added) {	// only use largest max_added overlaps
		std::sort(overlaps + sid, overlaps + eid, CompareOverlapByOverlapSize());
	}
	for (idx_t i(L); i < R; ++i) {
		Overlap& ovlp(overlaps[i]);
		qstr.resize(ovlp.qsize);
		reads.GetSequence(ovlp.qid, ovlp.qdir == FWD, qstr.data(), ovlp.qsize);
		const idx_t qext(ovlp.qdir == FWD ? ovlp.qext : ovlp.qsize - 1 - ovlp.qext);
		const idx_t sext(ovlp.sext);
		drd = drd_s;
		const bool r(GetAlignment(qstr.data(), qext, qstr.size(), tstr.data(), sext, tstr.size(), drd, m5, 0.15, min_align_size));
		if (r) {
			normalize_gaps(m5qaln(m5), m5saln(m5), strlen(m5qaln(m5)), nqstr, ntstr, true);
			meap_add_one_aln(nqstr, ntstr, m5soff(m5), cns_table, tstr.data());
			cns_vec.add_aln(m5soff(m5), m5send(m5), nqstr, ntstr);
		}
	}
	std::vector<MappingRange> mranges, eranges;
	cns_vec.get_mapping_ranges(mranges);
	get_effective_ranges(mranges, eranges, read_size, ctd.rco.min_size);
	consensus_worker(cns_table, pctd.id_list, cns_vec, nqstr, ntstr, eranges, ctd.rco.min_cov, ctd.rco.min_size, read_id, cns_results);
}

void
consensus_one_read_m4_nanopore(ConsensusThreadData& ctd, ConsensusPerThreadData &pctd, const idx_t read_id, const idx_t sid, const idx_t eid)
{
	PackedDB& reads = ctd.reads;
	ExtensionCandidate* overlaps = (ExtensionCandidate*)pctd.candidates;
	DiffRunningData* drd_s = pctd.drd_s;
	DiffRunningData* drd = NULL;
	M5Record& m5 = pctd.m5;
	CnsAlns& cns_vec = pctd.cns_alns;
	std::vector<CnsResult>& cns_results = pctd.cns_results;
	const idx_t read_size = overlaps[read_id].ssize;
	std::vector<char>& qstr = pctd.query;
	std::vector<char>& tstr = pctd.target;
	tstr.resize(read_size);
	reads.GetSequence(read_id, true, tstr.data(), read_size);
	std::string& nqstr = pctd.qaln;
	std::string& ntstr = pctd.saln;
	const int min_align_size = ctd.rco.min_align_size;
	const double min_mapping_ratio = ctd.rco.min_mapping_ratio - 0.02;

	idx_t L, R;
	if (eid - sid <= MAX_CNS_OVLPS)
	{
		L = sid;
		R = eid;
	}
	else
	{
		L = sid;
		R = L + MAX_CNS_OVLPS;    
		std::sort(overlaps + sid, overlaps + eid, CompareOverlapByOverlapSize());
	}

	CnsTableItem* cns_table = pctd.cns_table;
	std::for_each(cns_table, cns_table + read_size, CnsTableItemCleaner());
	cns_vec.clear();
	for (idx_t i = L; i < R; ++i)
	{
		Overlap& ovlp = overlaps[i];
		qstr.resize(ovlp.qsize);
		reads.GetSequence(ovlp.qid, ovlp.qdir == FWD, qstr.data(), ovlp.qsize);
		idx_t qext = ovlp.qext;
		idx_t sext = ovlp.sext;
		if (ovlp.qdir == REV) qext = ovlp.qsize - 1 - qext;
		drd = drd_s;
		bool r = GetAlignment(qstr.data(), qext, qstr.size(), tstr.data(), sext, tstr.size(), drd, m5, 0.20, min_align_size);
		if (r && check_ovlp_mapping_range(m5qoff(m5), m5qend(m5), ovlp.qsize, m5soff(m5), m5send(m5), ovlp.ssize, min_mapping_ratio))
		{
			normalize_gaps(m5qaln(m5), m5saln(m5), strlen(m5qaln(m5)), nqstr, ntstr, true);
			meap_add_one_aln(nqstr, ntstr, m5soff(m5), cns_table, tstr.data());
			cns_vec.add_aln(m5soff(m5), m5send(m5), nqstr, ntstr);
		}
	}
	
	std::vector<MappingRange> mranges, eranges;
	eranges.push_back(MappingRange(0, read_size));

	consensus_worker(cns_table, pctd.id_list, cns_vec, nqstr, ntstr, eranges, ctd.rco.min_cov, ctd.rco.min_size, read_id,  cns_results);
}

inline bool
check_cov_stats(u1_t* cov_stats, int soff, int send)
{
	const int max_cov = 20;
	int n = 0;
	for (int i = soff; i < send; ++i)
		if (cov_stats[i] >= max_cov)
			++n;
	if (send - soff >= n + 200) 
	{
		for (int i = soff; i < send; ++i) ++cov_stats[i];
		return true;
	}
	return false;
}

void consensus_one_read_can_pacbio(ConsensusThreadData& ctd, ConsensusPerThreadData& pctd, const idx_t read_id, const idx_t sid, idx_t eid) {
	const PackedDB& reads(ctd.reads);
	ExtensionCandidateCompressed* candidates((ExtensionCandidateCompressed*)pctd.candidates);
	DiffRunningData* const drd_s(pctd.drd_s);
	M5Record& m5(pctd.m5);
	CnsAlns& cns_vec(pctd.cns_alns);
	std::vector<CnsResult>& cns_results(pctd.cns_results);
	const idx_t ssize(reads.read_size(read_id));
	std::vector<char>& qstr(pctd.query);
	std::vector<char>& tstr(pctd.target);
	tstr.resize(ssize);
	reads.GetSequence(read_id, true, tstr.data(), ssize);
	std::string& nqstr(pctd.qaln);
	std::string& ntstr(pctd.saln);
	const int min_align_size(ctd.rco.min_align_size);
	const double min_mapping_ratio(ctd.rco.min_mapping_ratio - 0.02);
	int num_added(0);
	const int max_added(60);
	eid = std::min(eid, sid + 200);		// max of 200 extents
	CnsTableItem* cns_table(pctd.cns_table);
	std::for_each(cns_table, cns_table + ssize, CnsTableItemCleaner());
	cns_vec.clear();
	std::set<idx_t> used_ids;
	u1_t* cov_stats(pctd.id_list);
	std::fill(cov_stats, cov_stats + MAX_SEQ_SIZE, 0);
	for (idx_t i(sid); i < eid && num_added < max_added; ++i) {
		const ExtensionCandidateCompressed& ec(candidates[i]);
		if (used_ids.find(ec.qid) != used_ids.end()) {
			continue;
		}
		const idx_t qsize(reads.read_size(ec.qid));
		qstr.resize(qsize);
		reads.GetSequence(ec.qid, ec.qdir() == FWD, qstr.data(), qsize);
		const idx_t qext(ec.qdir() == FWD ? ec.qext() : qsize - 1 - ec.qext());
		const bool r(GetAlignment(qstr.data(), qext, qsize, tstr.data(), ec.sext, tstr.size(), drd_s, m5, 0.15, min_align_size));
		if (r && check_ovlp_mapping_range(m5qoff(m5), m5qend(m5), qsize, m5soff(m5), m5send(m5), ssize, min_mapping_ratio)) {
			if (check_cov_stats(cov_stats, m5soff(m5), m5send(m5))) {
				++num_added;
				used_ids.insert(ec.qid);
				normalize_gaps(m5qaln(m5), m5saln(m5), strlen(m5qaln(m5)), nqstr, ntstr, true);
				meap_add_one_aln(nqstr, ntstr, m5soff(m5), cns_table, tstr.data());
				cns_vec.add_aln(m5soff(m5), m5send(m5), nqstr, ntstr);
			}
		}
	}
	std::vector<MappingRange> mranges, eranges;
	cns_vec.get_mapping_ranges(mranges);
	get_effective_ranges(mranges, eranges, ssize, ctd.rco.min_size);
	if (ctd.rco.full_reads) {
		consensus_worker_one_read(cns_table, pctd.id_list, cns_vec, nqstr, ntstr, eranges, ctd.rco.min_cov, ctd.rco.min_size, read_id, tstr, cns_results);
	} else {
		consensus_worker(cns_table, pctd.id_list, cns_vec, nqstr, ntstr, eranges, ctd.rco.min_cov, ctd.rco.min_size, read_id, cns_results);
	}
}

void consensus_one_read_can_nanopore(ConsensusThreadData& ctd, ConsensusPerThreadData &pctd, const idx_t read_id, const idx_t sid, const idx_t eid) {
	PackedDB& reads(ctd.reads);
	ExtensionCandidateCompressed* candidates((ExtensionCandidateCompressed*)pctd.candidates);
	DiffRunningData* const drd_s(pctd.drd_s);
	DiffRunningData* drd(NULL);
	M5Record& m5(pctd.m5);
	CnsAlns& cns_vec(pctd.cns_alns);
	std::vector<CnsResult>& cns_results(pctd.cns_results);
	const idx_t ssize(reads.read_size(read_id));
	std::vector<char>& qstr = pctd.query;
	std::vector<char>& tstr = pctd.target;
	tstr.resize(ssize);
	reads.GetSequence(read_id, true, tstr.data(), ssize);
	std::string& nqstr(pctd.qaln);
	std::string& ntstr(pctd.saln);
	const int min_align_size(ctd.rco.min_align_size);
	const double min_mapping_ratio(ctd.rco.min_mapping_ratio - 0.02);
	int num_added(0);
	int num_ext(0);
	const int max_ext(200);
	CnsTableItem* cns_table(pctd.cns_table);
	std::for_each(cns_table, cns_table + ssize, CnsTableItemCleaner());
	cns_vec.clear();
	std::set<int> used_ids;
	u1_t* const cov_stats(pctd.id_list);
	std::fill(cov_stats, cov_stats + MAX_SEQ_SIZE, 0);
	for (idx_t i(sid); i < eid && num_added < MAX_CNS_OVLPS && num_ext < max_ext; ++i) {
		++num_ext;
		const ExtensionCandidateCompressed& ec(candidates[i]);
		if (used_ids.find(ec.qid) != used_ids.end()) {
			continue;
		}
		const idx_t qsize(reads.read_size(ec.qid));
		qstr.resize(qsize);
		reads.GetSequence(ec.qid, ec.qdir() == FWD, qstr.data(), qsize);
		const idx_t sext(ec.sext);
		const idx_t qext(ec.qdir() == FWD ? ec.qext() : qsize - 1 - ec.qext());
		drd = drd_s;
		const bool r(GetAlignment(qstr.data(), qext, qstr.size(), tstr.data(), sext, tstr.size(), drd, m5, 0.20, min_align_size));
		if (r && check_ovlp_mapping_range(m5qoff(m5), m5qend(m5), qsize, m5soff(m5), m5send(m5), ssize, min_mapping_ratio)) {
			if (check_cov_stats(cov_stats, m5soff(m5), m5send(m5))) {
				++num_added;
				used_ids.insert(ec.qid);
				normalize_gaps(m5qaln(m5), m5saln(m5), strlen(m5qaln(m5)), nqstr, ntstr, true);
				meap_add_one_aln(nqstr, ntstr, m5soff(m5), cns_table, tstr.data());
				cns_vec.add_aln(m5soff(m5), m5send(m5), nqstr, ntstr);
			}
		}
	}
	std::vector<MappingRange> mranges, eranges;
	eranges.push_back(MappingRange(0, ssize));
	consensus_worker(cns_table, pctd.id_list, cns_vec, nqstr, ntstr, eranges, ctd.rco.min_cov, ctd.rco.min_size, read_id,  cns_results);
}

} // namespace ns_meap_cns {
