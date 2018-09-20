#ifndef SARLACC_H
#define SARLACC_H

#include "Rcpp.h"

#include <stdexcept>
#include <vector>
#include <string>
#include <sstream>

extern "C" {

#include "Biostrings_interface.h"

}

extern "C" {

SEXP count_gaps_by_base(SEXP, SEXP, SEXP);
SEXP count_gaps_by_align(SEXP, SEXP, SEXP);

SEXP mask_bad_bases(SEXP, SEXP, SEXP, SEXP);
SEXP unmask_bases(SEXP, SEXP);

SEXP create_consensus_basic(SEXP, SEXP, SEXP);
SEXP create_consensus_basic_loop(SEXP, SEXP, SEXP);
SEXP create_consensus_quality(SEXP, SEXP, SEXP);
SEXP create_consensus_quality_loop(SEXP, SEXP, SEXP);

SEXP umi_group(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP compute_lev_masked(SEXP);

SEXP find_homopolymers(SEXP);
SEXP match_homopolymers(SEXP, SEXP);
SEXP find_errors(SEXP, SEXP);

SEXP get_aligned_sequence(SEXP, SEXP, SEXP, SEXP, SEXP);

}

#endif
