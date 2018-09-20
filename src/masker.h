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

class unmasker {
public:    
    unmasker(size_t);
    Rcpp::String unmask(const char*, size_t, const char*, size_t);
private:
    std::vector<char> buffer;
};


#endif
