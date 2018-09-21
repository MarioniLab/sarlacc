#include "sarlacc.h"
#include "DNA_input.h"
#include "masker.h"
#include "utils.h"

SEXP mask_bad_bases (SEXP sequences, SEXP qualities, SEXP encoding, SEXP threshold) {
    BEGIN_RCPP

    // Checking inputs.
    auto all_seq=process_DNA_input(sequences);
    const size_t nseq=all_seq->size();
    auto qholder=hold_XStringSet(qualities);

    if (nseq!=qholder.length) {
        throw std::runtime_error("sequence and quality vectors should have the same length");
    }
    
    masker seqmask(
        check_numeric_scalar(threshold, "quality threshold"),
        encoding
    );

    // Iterating through the sequences and masking bad bases.
    Rcpp::StringVector output(nseq);
    std::vector<char> buffer(10000, '\0');

    for (size_t i=0; i<nseq; ++i) {
        auto curpair=all_seq->get(i);
        const char* sstr=curpair.first;
        const size_t slen=curpair.second;

        auto curqual=get_elt_from_XStringSet_holder(&qholder, i);
        if (slen!=curqual.length) {
            throw std::runtime_error("sequence and quality strings should have the same length");
        }

        if (slen <= buffer.size()) {
            buffer.resize(slen+1);
        }

        seqmask.mask(slen, sstr, curqual.ptr, buffer.data());
        buffer[slen]='\0';
        output[i]=Rcpp::String(buffer.data());
    } 

    return output;
    END_RCPP
}
