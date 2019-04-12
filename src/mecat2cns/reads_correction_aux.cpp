#include "reads_correction_aux.h"

void normalize_gaps(const char* qstr, const char* tstr, const idx_t aln_size, std::string& qnorm, std::string& tnorm, const bool push)
{
    qnorm.clear();
    tnorm.clear();
    const char kGap = '-';

#ifndef NDEBUG
    int qcnt = 0, tcnt = 0;
    for (idx_t i = 0; i < aln_size; ++i)
    {
        const char qc = qstr[i];
        const char tc = tstr[i];
        if (qc != kGap) ++qcnt;
        if (tc != kGap) ++tcnt;
    }
#endif

    // convert mismatches to indels
    for (idx_t i = 0; i < aln_size; ++i)
    {
        const char qc = qstr[i];
        const char tc = tstr[i];
        if (qc != tc && qc != kGap && tc != kGap)
        { qnorm += kGap; qnorm += qc; tnorm += tc; tnorm += kGap; }
        else
        { qnorm += qc; tnorm += tc; }
    }

    // push gaps to the right, but not pass the end
    if (push)
    {
        idx_t qlen = qnorm.size();
        idx_t tlen = tnorm.size();
        for (idx_t i = 0; i < qlen - 1; ++i)
        {
            // push target gaps
            if (tnorm[i] == kGap)
            {
                idx_t j = i;
                while (1)
                {
                    const char c = tnorm[++j];
                    if (c != kGap || j > qlen - 1)
                    {
                        if (c == qnorm[i]) { tnorm[i] = c; tnorm[j] = kGap; }
                        break;
                    }
                }
            }
            // push query gaps
            if (qnorm[i] == kGap)
            {
                idx_t j = i;
                while (1)
                {
                    const char c = qnorm[++j];
                    if (c != kGap || j > tlen - 1)
                    {
                        if (c == tnorm[i]) { qnorm[i] = c; qnorm[j] = kGap; }
                        break;
                    }
                }
            }
        }
    }
    r_assert(qnorm.size() == tnorm.size());

#ifndef NDEBUG
    int qcnt2 = 0, tcnt2 = 0;
    for (std::string::const_iterator citer = qnorm.begin(); citer != qnorm.end(); ++citer)
        if ((*citer) != kGap) ++qcnt2;
    for (std::string::const_iterator citer = tnorm.begin(); citer != tnorm.end(); ++citer)
        if ((*citer) != kGap) ++tcnt2;
    d_assert(qcnt == qcnt2);
    d_assert(tcnt == tcnt2);
#endif
}

void allocate_ecs(ConsensusThreadData& data, ExtensionCandidate* const ec_list, const idx_t nec) {
	const int n(data.rco.num_threads);
	// split by number of ec's, rather than reads, since reads ids
	// are not contiguous and we could get empty lists
	for (idx_t i(0), k(0); k != n; ++k) {
		const idx_t start(i);
		// drop fractions here, as we'll likely add a few more ec's below
		i += (nec - i) / (n - k);
		if (i != nec) {			// include all ec's for the last read
			const int final_sid(ec_list[i].sid);
			for (++i; i != nec && ec_list[i].sid == final_sid; ++i) { }
		}
		data.data[k].num_candidates = i - start;
		data.data[k].candidates = ec_list + start;
	}
}

void allocate_ecs(ConsensusThreadData& data, ExtensionCandidateCompressed* const ec_list, const idx_t nec) {
	const int n(data.rco.num_threads);
	// split by number of ec's, rather than reads, since reads ids
	// are not contiguous and we could get empty lists
	for (idx_t i(0), k(0); k != n; ++k) {
		const idx_t start(i);
		// drop fractions here, as we'll likely add a few more ec's below
		i += (nec - i) / (n - k);
		if (i != nec) {			// include all ec's for the last read
			const int final_sid(ec_list[i].sid);
			for (++i; i != nec && ec_list[i].sid == final_sid; ++i) { }
		}
		data.data[k].num_candidates = i - start;
		data.data[k].candidates = ec_list + start;
	}
}
