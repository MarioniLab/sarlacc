#include "sarlacc.h"

#include "reference_align.h"
#include "utils.h"
#include "DNA_input.h"

#include <deque>
#include <stdexcept>

SEXP general_align(SEXP inputseq, SEXP inputqual, SEXP encoding, SEXP gapopen, SEXP gapext, SEXP reference, SEXP edit_only) {
    BEGIN_RCPP

    auto ref_seq=check_string(reference, "reference sequence");
    reference_align RA(ref_seq.size(), ref_seq.c_str(), encoding,
            check_numeric_scalar(gapopen, "gap opening penalty"), 
            check_numeric_scalar(gapext, "gap extension penalty")
    );

    auto sholder=process_DNA_input(inputseq);
    auto qholder=hold_XStringSet(inputqual);
    const size_t nseq=sholder->size();
    if (nseq!=qholder.length) {
        throw std::runtime_error("sequence and quality vectors should have the same length");
    }

    bool only_edit=check_logical_scalar(edit_only, "edit-only specification");

    // Setting up the output.
    Rcpp::NumericVector scores(nseq);
    Rcpp::IntegerVector edit_dist(nseq);
    Rcpp::StringVector refalign(only_edit ? 0 : nseq), qalign(only_edit ? 0 : nseq);
    std::vector<char> tmpref, tmpquery;

    for (size_t i=0; i<nseq; ++i) {
        auto curpair=sholder->get(i);
        const char* sstr=curpair.first;
        const size_t slen=curpair.second;

        auto curqual=get_elt_from_XStringSet_holder(&qholder, i);
        if (slen!=curqual.length) {
            throw std::runtime_error("sequence and quality strings should have the same length");
        }

        scores[i]=RA.align(slen, sstr, curqual.ptr, false);
        RA.fill_strings(tmpref, tmpquery, sstr);

        int& edist=edit_dist[i];
        for (size_t j=0; j<tmpref.size(); ++j) {
            if (tmpref[j]!=tmpquery[j]) {
                ++edist;
            }
        }

        if (!only_edit) {
            refalign[i]=tmpref.data();
            qalign[i]=tmpquery.data();
        }
    }

    return Rcpp::List::create(scores, edit_dist, refalign, qalign);
    END_RCPP
}
