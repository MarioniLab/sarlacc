#include "sarlacc.h"

/* The adjust_basepos_for_gaps function _ignores_ gaps, i.e., posstart and posend are defined on the bases, and are scalars.
 * The adjust_alignpos_for_gaps function includes gaps, i.e., posstart and posend are defined on the alignment string and are vectors.
 */

SEXP adjust_basepos_for_gaps(SEXP adapt_align, SEXP posstart, SEXP posend) {
    BEGIN_RCPP
    Rcpp::StringVector alignments(adapt_align);
    const size_t N=alignments.size();
    const size_t Start=check_integer_scalar(posstart, "UMI start position") - 1, 
                 End=check_integer_scalar(posend, "UMI end position"); // Zero indexed start, open end.

    Rcpp::IntegerVector beforeS(N), beforeE(N);
    for (size_t i=0; i<N; ++i) {
        Rcpp::String curstring(alignments[i]);
        const char* current=curstring.get_cstring();

        int counter=0;
        int& bumpStart=beforeS[i];
        int& bumpEnd=beforeE[i];
        
        while (*current!='\0' && counter < End) {
            if (*current!='-') {
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

SEXP adjust_alignpos_for_gaps(SEXP read_align, SEXP posstart, SEXP posend) {
    BEGIN_RCPP
    Rcpp::StringVector alignments(read_align);
    const size_t N=alignments.size();

    // Zero indexed start, open end.
    Rcpp::IntegerVector Starts(posstart), Ends(posend);
    if (N!=Starts.size() || N!=Ends.size()) {
        throw std::runtime_error("read alignments are not the same length as start/end vectors"); 
    }

    Rcpp::IntegerVector beforeS(N), beforeE(N);
    for (size_t i=0; i<N; ++i) {
        Rcpp::String curstring(alignments[i]);
        const char* current=curstring.get_cstring();

        int counter=0;
        int& bumpStart=beforeS[i];
        int& bumpEnd=beforeE[i];

        int Start=Starts[i]-1; // zero indexed start, open end.
        int End=Ends[i];
        
        while (*current!='\0' && counter < End) {
            if (*current=='-') {
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

