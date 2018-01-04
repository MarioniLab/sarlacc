#include "Rcpp.h"
#include <stdexcept>

extern "C" {

SEXP mask_bad_bases(SEXP, SEXP, SEXP, SEXP);

SEXP count_deletions(SEXP, SEXP, SEXP);

}
