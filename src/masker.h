#ifndef MASKER_H
#define MASKER_H

#include "DNA_input.h"
#include "Rcpp.h"

class masker {
public:
    masker(double, size_t);
    Rcpp::String mask(const char*, size_t, double*);
private:
    std::vector<char> buffer;
    double threshold;
};

class unmasker {
public:    
    unmasker(size_t);
    Rcpp::String unmask(const char*, size_t, const char*, size_t);
private:
    std::vector<char> buffer;
};


#endif
