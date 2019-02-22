#ifndef REFERENCE_ALIGN_H
#define REFERENCE_ALIGN_H

#include "Rcpp.h"
#include <vector>
#include <deque>
#include "quality_encoding.h"

class reference_align {
public:
    reference_align(size_t, const char*, Rcpp::NumericVector, double, double);
    double align(size_t, const char*, const char*, bool=true);

    void backtrack(std::deque<size_t>&, std::deque<size_t>&) const;
    void backtrack(std::deque<char>&, std::deque<char>&, const char*) const;

    static std::pair<size_t, size_t> map(const std::deque<size_t>&, const std::deque<size_t>&, size_t, size_t); 
private:
    // Reference related.
    size_t rlen;
    const char* rseq;

    // Gap related.
    double gap_open, gap_ext;

    // Quality related.
    std::deque<std::vector<double> > precomputed_match, precomputed_mismatch;
    char offset;
    size_t available;
    void create_qualities(Rcpp::NumericVector);

    // DP matrix-related.
    bool aligned=false;
    size_t nrows=1000;
    std::deque<double> scores, affine_left;
    enum DIRECTION { up, left, diag };
    std::deque<DIRECTION> directions;

    void align_column(std::deque<DIRECTION>::iterator, char, size_t, const char*, const char*, bool);
    double compute_cost (char, char, char) const;
    double precomputed_cost (int, bool, char) const;
};

#endif
