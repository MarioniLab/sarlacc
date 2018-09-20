#include "sarlacc.h"
#include "DNA_input.h"
#include "masker.h"
#include "utils.h"

SEXP mask_bad_bases (SEXP sequences, SEXP qualities, SEXP encoding, SEXP threshold) {
    BEGIN_RCPP

    // Checking inputs.
    auto all_seq=process_DNA_input(sequences);
    const size_t nseq=all_seq->size();

    Rcpp::StringVector all_qual(qualities);
    if (nseq!=all_qual.size()) {
        throw std::runtime_error("sequence and quality vectors should have the same length");
    }
    
    masker seqmask(
        check_numeric_scalar(threshold, "quality threshold"),
        encoding
    );

    // Iterating through the sequences and masking bad bases.
    Rcpp::StringVector output(nseq);
    std::vector<char> buffer(100, '\0');

    for (size_t i=0; i<nseq; ++i) {
        auto curpair=all_seq->get(i);
        const char* sstr=curpair.first;
        const size_t slen=curpair.second;

        if (slen <= buffer.size()) {
            buffer.resize(slen+1);
        }

        Rcpp::String qualstr=all_qual[i];
        const char* qptr=qualstr.get_cstring();

        seqmask.mask(slen, sstr, qptr, buffer.data());
        buffer[slen]='\0';
        output[i]=Rcpp::String(buffer.data());
    } 

    return output;
    END_RCPP
}

/* This performs the reverse operation in order to obtain the 
 * alignment strings but with the low-quality bases unmasked.
 */

SEXP unmask_bases (SEXP alignments, SEXP originals) {
    BEGIN_RCPP

    auto all_aln=process_DNA_input(alignments);
    const size_t nseq=all_aln->size();
    auto all_seq=process_DNA_input(originals); 
    if (nseq!=all_seq->size()) {
        throw std::runtime_error("alignment and original sequences should have the same number of entries");
    }

    const size_t aln_width=check_alignment_width(all_aln.get());
  
    // Running through the output and restoring the masked bases. 
    Rcpp::StringVector output(nseq);
    unmasker sequnmask(aln_width);

    for (size_t i=0; i<nseq; ++i) {
        auto apair=all_aln->get(i);
        auto spair=all_seq->get(i);
        output[i]=sequnmask.unmask(apair.first, aln_width, spair.first, spair.second);
    }
    return output;
    END_RCPP
}
