#include "sarlacc.h"

#include "DNA_input.h"

#include <stdexcept>
#include <vector>

/* Obtain alignment strings but with the low-quality bases unmasked.
 * Requires the original sequences, and only works on global alignments.
 */

SEXP unmask_alignment(SEXP alignments, SEXP originals) {
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
    std::vector<char> buffer(aln_width+1, '\0');

    for (size_t i=0; i<nseq; ++i) {
        auto apair=all_aln->get(i);
        const char* masked=apair.first;

        auto spair=all_seq->get(i);
        const char * origin=spair.first;
        const size_t orilen=spair.second;

        size_t pos_nominal=0;
        for (size_t pos=0; pos<aln_width; ++pos) {
            char& outbase=(buffer[pos]=masked[pos]);

            if (outbase!='-') {
                if (outbase=='N' || outbase=='n') {
                    if (pos_nominal >= orilen) {
                        throw std::runtime_error("sequence in alignment string is longer than the original");
                    }
                    outbase=origin[pos_nominal];
                }
                ++pos_nominal;
            }
        }

        if (pos_nominal!=orilen) {
            throw std::runtime_error("original sequence and that in the alignment string have different lengths");
        }

        output[i]=Rcpp::String(buffer.data());
    }
    return output;
    END_RCPP
}
