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
    void backtrack(bool=false);
    std::pair<size_t, size_t> map(size_t, size_t) const; 
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
    enum DIRECTION { up, left, diag };
    typedef std::pair<DIRECTION, double> dpentry;
    size_t nrows=1000;
    bool aligned=false;
    std::deque<dpentry> dpmatrix;
    std::deque<double> affine_left;

    void align_column(std::deque<dpentry>::iterator, std::deque<double>::iterator, char, size_t, const char*, const char*, bool);
    double compute_cost (char, char, char) const;
    double precomputed_cost (int, bool, char) const;

    // Backtrack-related.
    bool backtracked=false;
    std::deque<size_t> backtrack_start, backtrack_end;
};

#endif
