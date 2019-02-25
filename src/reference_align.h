#ifndef REFERENCE_ALIGN_H
#define REFERENCE_ALIGN_H

#include "Rcpp.h"

#include "quality_encoding.h"

#include <vector>
#include <deque>

class reference_align {
public:
    reference_align(size_t, const char*, Rcpp::NumericVector, double, double);
    double align(size_t, const char*, const char*, bool=true);

    // Class to map each reference position to a range on the query.
    struct querymap {
        std::deque<std::pair<bool, size_t> > mapping;
        size_t nrows;
        std::pair<size_t, size_t> operator()(size_t, size_t, bool=false) const; 
    };
    void fill_map(querymap&) const;

    void fill_strings(std::vector<char>&, std::vector<char>&, const char*) const;
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
    std::deque<double> scores, left_jump_scores;
    std::deque<size_t> left_jump_points;
    std::deque<int> directions;

    void align_column(std::deque<int>::iterator, size_t, size_t, const char*, const char*, bool);
    double compute_cost (char, char, char) const;
    double precomputed_cost (int, bool, char) const;

    // Backtrack related.
    template <class U>
    void backtrack(U&) const;
};

#endif
