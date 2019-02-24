#ifndef SARLACC_H
#define SARLACC_H

#include "Rcpp.h"

extern "C" {

#include "Biostrings_interface.h"

}

extern "C" {

SEXP adaptor_align(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP adaptor_align_score_only(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP barcode_align(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP general_align(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
        
SEXP mask_bad_bases(SEXP, SEXP, SEXP, SEXP);
SEXP unmask_alignment(SEXP, SEXP);

SEXP create_consensus_basic(SEXP, SEXP, SEXP);
SEXP create_consensus_basic_loop(SEXP, SEXP, SEXP);
SEXP create_consensus_quality(SEXP, SEXP, SEXP, SEXP);
SEXP create_consensus_quality_loop(SEXP, SEXP, SEXP, SEXP);

SEXP umi_group(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP fast_levdist_test(SEXP, SEXP, SEXP);
SEXP cluster_umis_test(SEXP);
SEXP compute_lev_masked(SEXP);

SEXP find_homopolymers(SEXP);
SEXP match_homopolymers(SEXP, SEXP);
SEXP find_errors(SEXP, SEXP);

SEXP quick_msa(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

}

#endif
