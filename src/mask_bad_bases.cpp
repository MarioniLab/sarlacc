#include "sarlacc.h"

#include "quality_encoding.h"
#include "utils.h"
#include "DNA_input.h"

#include <stdexcept>
#include <vector>

SEXP mask_bad_bases (SEXP sequences, SEXP qualities, SEXP encoding, SEXP threshold) {
    BEGIN_RCPP

    // Checking inputs.
    auto sholder=process_DNA_input(sequences);
    const size_t nseq=sholder->size();
    auto qholder=hold_XStringSet(qualities);

    if (nseq!=qholder.length) {
        throw std::runtime_error("sequence and quality vectors should have the same length");
    }
    
    quality_encoding seqmask(encoding);
    double maxerr=check_numeric_scalar(threshold, "quality threshold");

    // Iterating through the sequences and masking bad bases.
    Rcpp::StringVector output(nseq);
    std::vector<char> buffer(10000, '\0');

    for (size_t i=0; i<nseq; ++i) {
        auto curpair=sholder->get(i);
        const char* sstr=curpair.first;
        const size_t slen=curpair.second;

        auto curqual=get_elt_from_XStringSet_holder(&qholder, i);
        if (slen!=curqual.length) {
            throw std::runtime_error("sequence and quality strings should have the same length");
        }

        if (slen <= buffer.size()) {
            buffer.resize(slen+1);
        }

        for (size_t counter=0; counter<slen; ++counter) {
            buffer[counter]=(seqmask.to_error(curqual.ptr[counter]) > maxerr ? 'N' : sstr[counter]);
        }
        buffer[slen]='\0';
        output[i]=Rcpp::String(buffer.data());
    } 

    return output;
    END_RCPP
}
