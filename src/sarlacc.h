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

SEXP mask_bad_bases(SEXP, SEXP, SEXP);

SEXP count_deletions(SEXP, SEXP, SEXP);

}
