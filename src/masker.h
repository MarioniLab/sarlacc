#ifndef MASKER_H
#define MASKER_H

#include "DNA_input.h"
#include "Rcpp.h"

class masker {
public:
    masker(double, Rcpp::NumericVector);
    void mask(size_t, const char*, const char*, char*);
private:
    char offset;
};

#endif
