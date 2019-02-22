#include "sarlacc.h"

#include "reference_align.h"
#include "utils.h"
#include "DNA_input.h"

#include <deque>
#include <stdexcept>

SEXP barcode_align(SEXP barcodeseq, SEXP barcodequal, SEXP encoding, SEXP gapopen, SEXP gapext, SEXP reference) {
    BEGIN_RCPP

    auto ref_barcode_seq=check_string(reference, "barcode sequence");
    reference_align RA(ref_barcode_seq.size(), ref_barcode_seq.c_str(), encoding,
            check_numeric_scalar(gapopen, "gap opening penalty"), 
            check_numeric_scalar(gapext, "gap extension penalty")
    );

    auto sholder=process_DNA_input(barcodeseq);
    auto qholder=hold_XStringSet(barcodequal);
    const size_t nseq=sholder->size();
    if (nseq!=qholder.length) {
        throw std::runtime_error("sequence and quality vectors should have the same length");
    }

    // Setting up the output.
    Rcpp::NumericVector scores(nseq);

    for (size_t i=0; i<nseq; ++i) {
        auto curpair=sholder->get(i);
        const char* sstr=curpair.first;
        const size_t slen=curpair.second;
        
        auto curqual=get_elt_from_XStringSet_holder(&qholder, i);
        if (slen!=curqual.length) {
            throw std::runtime_error("sequence and quality strings should have the same length");
        }

        scores[i]=RA.align(slen, sstr, curqual.ptr, false);
    }

    return scores;
    END_RCPP
}
