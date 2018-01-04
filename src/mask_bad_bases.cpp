#include "sarlacc.h"

SEXP mask_bad_bases (SEXP sequences, SEXP qualities, SEXP threshold) {
    BEGIN_RCPP

    // Checking inputs.       
    auto seq=hold_XStringSet(sequences); 
    const size_t nseq=get_length_from_XStringSet_holder(&seq);
    
    auto qual=hold_XStringSet(qualities); 
    if (nseq!=get_length_from_XStringSet_holder(&qual)) { 
        throw std::runtime_error("sequence and quality vectors should have the same length");
    }

    std::string curstring=check_string(threshold, "quality threshold");
    if (curstring.size()!=1) {
        throw std::runtime_error("quality threshold should be a string of length 1");
    }
    const char lowerbound=curstring[0]; 
    
    // Getting the width of the buffer that should be set.
    int buffersize=0;
    {
        for (size_t i=0; i<nseq; ++i) { 
            auto curseq=get_elt_from_XStringSet_holder(&seq, i);
            if (curseq.length > buffersize) {
                buffersize=curseq.length;
            }
        }
        ++buffersize; // for the NULL.
    }

    // Iterating through the sequences and masking bad bases.
    Rcpp::StringVector output(nseq);
    std::vector<char> buffer(buffersize);
    
    for (size_t i=0; i<nseq; ++i) {
        auto curseq=get_elt_from_XStringSet_holder(&seq, i);
        const char* sstr=curseq.ptr;
        const size_t slen=curseq.length;

        auto curqual=get_elt_from_XStringSet_holder(&qual, i);
        const char* qstr=curqual.ptr;
        const size_t qlen=curqual.length;
        
        if (slen!=qlen) {
            throw std::runtime_error("sequence and quality strings are not the same length");
        }

        for (size_t counter=0; counter<slen; ++counter) {
            buffer[counter]=(*qstr < lowerbound ? 'N' : DNAdecode(*sstr));
            ++qstr;
            ++sstr;
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

    auto aln=hold_XStringSet(alignments); 
    const size_t nseq=get_length_from_XStringSet_holder(&aln);
    
    auto seq=hold_XStringSet(originals); 
    if (nseq!=get_length_from_XStringSet_holder(&seq)) { 
        throw std::runtime_error("alignment and original sequences should have the same length");
    }

    // Getting the width of the buffer that should be set.
    int buffersize=0;
    {
        for (size_t i=0; i<nseq; ++i) { 
            auto curaln=get_elt_from_XStringSet_holder(&aln, i);
            if (curaln.length > buffersize) {
                buffersize=curaln.length;
            }
        }
        ++buffersize; // for the NULL.
    }
  
    // Running through the output and restoring the masked bases. 
    Rcpp::StringVector output(nseq);
    std::vector<char> buffer(buffersize);
    
    for (size_t i=0; i<nseq; ++i) {
        auto curaln=get_elt_from_XStringSet_holder(&aln, i);
        const char* astr=curaln.ptr;
        const size_t alen=curaln.length;

        auto curseq=get_elt_from_XStringSet_holder(&seq, i);
        const char* sstr=curseq.ptr;
        const size_t slen=curseq.length;

        size_t a_nominal=0;
        for (size_t a=0; a<alen; ++a) {
            char& outbase=(buffer[a]=DNAdecode(astr[a]));
    
            if (outbase!='-') { 
                if (outbase=='N' || outbase=='n') {
                    if (a_nominal >= slen) {
                        throw std::runtime_error("sequence in alignment string is longer than the original");
                    }
                    outbase=DNAdecode(sstr[a_nominal]);
                }
                ++a_nominal;
            }
        }

        if (a_nominal!=slen) {
            throw std::runtime_error("original sequence and that in the alignment string have different lengths");
              
        }

        buffer[alen]='\0';
        output[i]=Rcpp::String(buffer.data());
    }
    return output;
    END_RCPP
}
