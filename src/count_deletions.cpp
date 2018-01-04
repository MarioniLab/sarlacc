#include "sarlacc.h"

SEXP count_deletions(SEXP _alignments, SEXP posstart, SEXP posend) {
    BEGIN_RCPP
    Rcpp::StringVector alignments(_alignments);
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
                    if (counter <= Start) {
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

