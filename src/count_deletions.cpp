#include "sarlacc.h"

SEXP count_deletions(SEXP _alignments, SEXP _posstart, SEXP _posend) {
    BEGIN_RCPP
    Rcpp::StringVector alignments(_alignments);
    const size_t N=alignments.size();
    
    Rcpp::IntegerVector posstart(_posstart);
    if (posstart.size()!=1) {
        throw std::runtime_error("UMI position start should be an integer scalar");
    }
    Rcpp::IntegerVector posend(_posend);
    if (posend.size()!=1) {
        throw std::runtime_error("UMI position end should be an integer scalar");
    }
    const size_t Start=posstart[0] - 1, End=posend[0]; // Zero indexed start, open end.

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

