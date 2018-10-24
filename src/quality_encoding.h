#ifndef QUALITY_ENCODING_H
#define QUALITY_ENCODING_H

#include "Rcpp.h"

class quality_encoding {
public:
    quality_encoding(Rcpp::NumericVector);
    char lowest() const;
    Rcpp::NumericVector get_errors() const;
    double to_error(char) const;
private:
    Rcpp::NumericVector errors;
    char offset;
};

#endif
