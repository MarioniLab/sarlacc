#include "sarlacc.h"

int check_integer_scalar(Rcpp::RObject, const char*);

double check_numeric_scalar(Rcpp::RObject, const char*);

bool check_logical_scalar(Rcpp::RObject, const char*);

std::string check_string(Rcpp::RObject, const char*);
