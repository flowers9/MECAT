#ifndef MEAP_CORRECTION_H
#define MEAP_CORRECTION_H

#include "reads_correction_aux.h"	// ConsensusPerThreadData, ConsensusThreadData
#include "../common/defs.h"	// idx_t

void consensus_one_read_m4_pacbio(ConsensusThreadData& ctd, ConsensusPerThreadData &cptd, const idx_t read_id, const idx_t sid, const idx_t eid);

void consensus_one_read_m4_nanopore(ConsensusThreadData& ctd, ConsensusPerThreadData &cptd, const idx_t read_id, const idx_t sid, const idx_t eid);

void consensus_one_read_can_pacbio(ConsensusThreadData& ctd, ConsensusPerThreadData &cptd, const idx_t read_id, const idx_t sid, const idx_t eid);

void consensus_one_read_can_nanopore(ConsensusThreadData& ctd, ConsensusPerThreadData &cptd, const idx_t read_id, const idx_t sid, const idx_t eid);

#endif // MEAP_CORRECTION_H
