#include "sarlacc.h"

SEXP mask_bad_bases (SEXP sequences, SEXP qualities, SEXP maxwidth, SEXP threshold) {
    BEGIN_RCPP

    // Checking inputs.        
    Rcpp::StringVector seq(sequences), qual(qualities);
    const size_t nseq=seq.size();
    if (nseq!=qual.size()) {
        throw std::runtime_error("sequence and quality vectors should have the same length");
    }

    char lowerbound=0;
    {
        Rcpp::StringVector limit(threshold);   
        if (limit.size()!=1) {
            throw std::runtime_error("threshold should be a single character");
        }
        Rcpp::String first(limit[0]);
        lowerbound=first.get_cstring()[0]; // should have at least 1 character, the NULL.
    }

    int buffersize=0;
    {
        Rcpp::IntegerVector mw(maxwidth);
        if (mw.size()!=1) {
            throw std::runtime_error("maximum string width should be an integer scalar");
        }
        buffersize=mw[0]+1; // for the NULL.
    }

    // Iterating through the sequences and masking bad bases.
    Rcpp::StringVector output(nseq);
    std::vector<char> buffer(buffersize);
    
    for (size_t i=0; i<nseq; ++i) {
        Rcpp::String curseq(seq[i]);
        Rcpp::String curqual(qual[i]);
        const char* sstr=curseq.get_cstring();
        const char* qstr=curqual.get_cstring();

        int counter=0;
        while (*sstr!='\0' && *qstr!='\0') {
            buffer[counter]=(*qstr < lowerbound ? 'N' : *sstr);
            ++qstr;
            ++sstr;
            ++counter;
            if (counter >= buffersize) {
                throw std::runtime_error("maximum string length exceeded");
            }
        }
        if (*sstr!=*qstr) {
            throw std::runtime_error("sequence and quality strings are not the same length");
        }

        output[i]=Rcpp::String(buffer.data());
    } 

    return output;
    END_RCPP
}
