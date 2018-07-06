#include "sarlacc.h"
#include "DNA_input.h"
#include "utils.h"

SEXP mask_bad_bases (SEXP sequences, SEXP qualities, SEXP threshold) {
    BEGIN_RCPP

    // Checking inputs.
    auto all_seq=process_DNA_input(sequences);
    const size_t nseq=all_seq->size();

    Rcpp::List all_qual(qualities);
    if (nseq!=all_qual.size()) {
        throw std::runtime_error("sequence and quality vectors should have the same length");
    }
    const double max_error=check_numeric_scalar(threshold, "quality threshold");
    
    // Getting the length of the buffer that should be set.
    int buffersize=get_max_width(all_seq.get());
    ++buffersize; // for the NULL.

    // Iterating through the sequences and masking bad bases.
    Rcpp::StringVector output(nseq);
    std::vector<char> buffer(buffersize);
    
    for (size_t i=0; i<nseq; ++i) {
        auto curpair=all_seq->get(i);
        const char* sstr=curpair.first;
        const size_t slen=curpair.second;

        Rcpp::NumericVector curqual(all_qual[i]);
        if (slen!=curqual.size()) {
            Rprintf("%s %i %i\n", sstr, slen, curqual.size());
            throw std::runtime_error("sequence and quality strings are not the same length");
        }

        for (size_t counter=0; counter<slen; ++counter, ++sstr) {
            buffer[counter]=(curqual[counter] > max_error ? 'N' : *sstr);
        }

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

    // Getting the width of the buffer that should be set (+1 for the NULL).
    const size_t aln_width=check_alignment_width(all_aln.get());
    const size_t buffersize=aln_width + 1;
  
    // Running through the output and restoring the masked bases. 
    Rcpp::StringVector output(nseq);
    std::vector<char> buffer(buffersize);
    
    for (size_t i=0; i<nseq; ++i) {
        auto apair=all_aln->get(i);
        const char* astr=apair.first;

        auto spair=all_seq->get(i);
        const char* sstr=spair.first;
        const size_t slen=spair.second;

        size_t a_nominal=0;
        for (size_t a=0; a<aln_width; ++a) {
            char& outbase=(buffer[a]=astr[a]);
    
            if (outbase!='-') { 
                if (outbase=='N' || outbase=='n') {
                    if (a_nominal >= slen) {
                        throw std::runtime_error("sequence in alignment string is longer than the original");
                    }
                    outbase=sstr[a_nominal];
                }
                ++a_nominal;
            }
        }

        if (a_nominal!=slen) {
            throw std::runtime_error("original sequence and that in the alignment string have different lengths");
              
        }

        buffer[aln_width]='\0';
        output[i]=Rcpp::String(buffer.data());
    }
    return output;
    END_RCPP
}
