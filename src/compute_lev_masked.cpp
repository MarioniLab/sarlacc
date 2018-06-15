#include "sarlacc.h"

/* This function computes a masked version of the Levenshtein distance,
 * where N's always contribute a mismatch no matter what they are matched to.
 * In effect, they are treated as 'missing' values.
 */

SEXP compute_lev_masked(SEXP seqs) {
    BEGIN_RCPP

    auto seq=hold_XStringSet(seqs);
    const size_t nseq=get_length_from_XStringSet_holder(&seq);
    Rcpp::NumericVector output((nseq-1)*nseq/2);
    auto oIt=output.begin();

    // Setting up the odds and ends required
    std::vector<double> col, prev_col;
    col.reserve(50);
    prev_col.reserve(50);

    for (size_t i=0; i<nseq; ++i) {
        auto iseq=get_elt_from_XStringSet_holder(&seq, i);
        const char* istr=iseq.ptr;
        const size_t ilen=iseq.length;

        if (prev_col.size() < ilen) {
            col.resize(ilen + 1);
            prev_col.resize(ilen + 1);
        }

        for (size_t j=i+1; j<nseq; ++j) {
            auto jseq=get_elt_from_XStringSet_holder(&seq, j);
            const char* jstr=jseq.ptr;
            const size_t jlen=jseq.length;

            // Computing the levenshtein distance.
            std::iota(prev_col.begin(), prev_col.end(), 0);
            for (size_t jx=0; jx<jlen; ++jx) {
                col[0]=jx+1;
                const char jbase=DNAdecode(jstr[jx]);

                for (size_t ix=0; ix<ilen; ++ix) {
                    const char ibase=DNAdecode(istr[ix]);
                    const double match_score = (jbase=='N' || ibase=='N' ? 0.5 : (jbase==ibase ? 0 : 1));
                    col[ix+1] = std::min({ prev_col[ix+1] + 1, col[ix] + 1, prev_col[ix] + match_score });
                }
                col.swap(prev_col);
            }

            (*oIt)=prev_col[ilen];
            ++oIt;
        }
    }

    return output;
    END_RCPP
}
