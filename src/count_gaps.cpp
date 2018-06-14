#include "sarlacc.h"
#include "utils.h"
#include "DNA_input.h"

/* The count_gaps_by_base function ignores gaps in the position specification, i.e., posstart and posend refer to base positions (and are scalars).
 * The count_gaps_by_align function includes gaps in the position specification, i.e., posstart and posend refer to alignment columns (and are vectors).
 */

SEXP count_gaps_by_base(SEXP adapt_align, SEXP posstart, SEXP posend) {
    BEGIN_RCPP
    auto alignments=process_DNA_input(adapt_align);
    const size_t N=alignments->size();

    const int Start=check_integer_scalar(posstart, "UMI start position") - 1, 
              End=check_integer_scalar(posend, "UMI end position"); // Zero indexed start, open end.
    if (Start<0 || End<0) {
        throw std::runtime_error("positions must be positive integers"); // message for R-level users, hence not 'non-negative'.
    }
    if (Start > End) {
        throw std::runtime_error("start position must be less than or equal to end position");
    }

    Rcpp::IntegerVector beforeS(N), beforeE(N);
    for (size_t i=0; i<N; ++i) {
        alignments->choose(i);
        const size_t len=alignments->length();
        const char* current=alignments->cstring();
        
        int counter=0;
        int& bumpStart=beforeS[i];
        int& bumpEnd=beforeE[i];
        
        while (counter < End) {
            if (counter > len) {
                throw std::runtime_error("end position exceeds the number of bases in the alignment string");
            }
            if (alignments->decode(*current)!='-') {
                ++counter;
            } else {
                if (counter < End) {
                    ++bumpEnd;
                    if (counter <= Start) { // '<=' ensures gaps just before the 'Start' base are included in 'bumpStart'.
                        ++bumpStart;
                    }
                }
            }
            ++current;
        } 
    }

    return Rcpp::List::create(beforeS, beforeE);
    END_RCPP
}

SEXP count_gaps_by_align(SEXP read_align, SEXP posstart, SEXP posend) {
    BEGIN_RCPP
    auto alignments=process_DNA_input(read_align);
    const size_t N=alignments->size();

    // Zero indexed start, open end.
    Rcpp::IntegerVector Starts(posstart), Ends(posend);
    if (N!=Starts.size() || N!=Ends.size()) {
        throw std::runtime_error("read alignments are not the same length as start/end vectors"); 
    }

    Rcpp::IntegerVector beforeS(N), beforeE(N);
    for (size_t i=0; i<N; ++i) {
        alignments->choose(i);
        const char* current=alignments->cstring();
        const size_t len=alignments->length();

        int Start=Starts[i]-1; // zero indexed start, open end.
        int End=Ends[i];
        if (Start<0 || End<0) {
            throw std::runtime_error("positions must be positive integers"); // message for R-level users, hence not 'non-negative'.
        }
        if (Start > End) {
            throw std::runtime_error("start position must be less than or equal to end position");
        }
        if (End > len) {
            throw std::runtime_error("end position exceeds the number of positions in the alignment string");
        }

        int counter=0;
        int& bumpStart=beforeS[i];
        int& bumpEnd=beforeE[i];

        while (counter < End) {
            if (alignments->decode(*current)=='-') {
                if (counter < End) {
                    ++bumpEnd;
                    if (counter < Start) { // NOT '<=', to avoid counting gaps at the 'Start' alignment position.
                        ++bumpStart;
                    }
                }
            }
            ++counter;
            ++current;
        } 
    }

    return Rcpp::List::create(beforeS, beforeE);
    END_RCPP
}

