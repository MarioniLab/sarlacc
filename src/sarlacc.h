#include "Rcpp.h"

#include <stdexcept>
#include <vector>
#include <string>
#include <sstream>

extern "C" {

#include "Biostrings_interface.h"

}

int check_integer_scalar(Rcpp::RObject, const char*);

double check_numeric_scalar(Rcpp::RObject, const char*);

bool check_logical_scalar(Rcpp::RObject, const char*);

std::string check_string(Rcpp::RObject, const char*);

extern "C" {

SEXP count_gaps_by_base(SEXP, SEXP, SEXP);
SEXP count_gaps_by_align(SEXP, SEXP, SEXP);

SEXP mask_bad_bases(SEXP, SEXP, SEXP);
SEXP unmask_bases(SEXP, SEXP);

SEXP create_consensus_basic(SEXP, SEXP, SEXP);
SEXP create_consensus_basic_loop(SEXP, SEXP, SEXP);
SEXP create_consensus_quality(SEXP, SEXP, SEXP);
SEXP create_consensus_quality_loop(SEXP, SEXP, SEXP);

SEXP umi_group(SEXP, SEXP, SEXP);
SEXP compute_lev_masked(SEXP);

SEXP get_kmer_matrix(SEXP, SEXP);

SEXP find_homopolymers(SEXP);
SEXP match_homopolymers(SEXP, SEXP);
SEXP find_errors(SEXP, SEXP);

SEXP get_aligned_sequence(SEXP, SEXP, SEXP, SEXP, SEXP);

}
