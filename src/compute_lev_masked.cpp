#include "sarlacc.h"

#include "DNA_input.h"

#include <vector>
#include <algorithm>

/* This function computes a masked version of the Levenshtein distance,
 * where N's always contribute a mismatch no matter what they are matched to.
 * In effect, they are treated as 'missing' values.
 */

SEXP compute_lev_masked(SEXP sequences) {
    BEGIN_RCPP

    auto seqs=process_DNA_input(sequences);
    const size_t nseq=seqs->size();
    Rcpp::NumericVector output((nseq-1)*nseq/2);
    auto oIt=output.begin();

    // Setting up the odds and ends required
    std::vector<double> col, prev_col;
    col.reserve(50);
    prev_col.reserve(50);

    for (size_t i=0; i<nseq; ++i) {
        auto iseq=seqs->get(i);
        const char* istr=iseq.first;
        const size_t ilen=iseq.second;

        if (prev_col.size() < ilen) {
            col.resize(ilen + 1);
            prev_col.resize(ilen + 1);
        }

        auto altseqs=process_DNA_input(sequences);

        for (size_t j=i+1; j<nseq; ++j) {
            auto jseq=altseqs->get(j);
            const char* jstr=jseq.first;
            const size_t jlen=jseq.second;

            // Computing the levenshtein distance.
            std::iota(prev_col.begin(), prev_col.end(), 0);
            for (size_t jx=0; jx<jlen; ++jx) {
                col[0]=jx+1;
                const char jbase=jstr[jx];

                for (size_t ix=0; ix<ilen; ++ix) {
                    const char ibase=istr[ix];
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
