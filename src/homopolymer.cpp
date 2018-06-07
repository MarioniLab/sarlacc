#include "sarlacc.h"

SEXP find_homopolymers (SEXP sequences) {
    BEGIN_RCPP

    // Checking inputs.       
    auto seq=hold_XStringSet(sequences); 
    const size_t nseq=get_length_from_XStringSet_holder(&seq);
    
    // Creating outputs.
    std::deque<int> collected_index, collected_pos, collected_size;
    std::deque<char> collected_base;
    
    // Running through each string and returning the number of times we see a homopolymer.
    for (size_t i=0; i<nseq; ++i) {
        auto curseq=get_elt_from_XStringSet_holder(&seq, i);
        const char* sstr=curseq.ptr;
        const size_t slen=curseq.length;

        if (slen) {
            size_t last_pos=0, cur_pos=0;
            char last_base=DNAdecode(sstr[0]);
            while ((++cur_pos) < slen) {
                char outbase=DNAdecode(sstr[cur_pos]);

                if (outbase==last_base) {
                    while ((++cur_pos) < slen) {
                        outbase=DNAdecode(sstr[cur_pos]);
                        if (outbase!=last_base) { 
                            break;
                        }
                    }

                    collected_index.push_back(i);
                    collected_pos.push_back(last_pos+1);
                    collected_size.push_back(cur_pos-last_pos);
                    collected_base.push_back(last_base);
                }

                // This is always the first element where out_base!=last_base.
                last_pos=cur_pos;
                last_base=outbase;
            }
        }
    }

    Rcpp::List output(4);
    output[0]=Rcpp::IntegerVector(collected_index.begin(), collected_index.end());
    output[1]=Rcpp::IntegerVector(collected_pos.begin(), collected_pos.end());
    output[2]=Rcpp::IntegerVector(collected_size.begin(), collected_size.end());

    // Bit more effort to convert character strings.
    Rcpp::StringVector bases(collected_base.size());
    auto bIt=bases.begin();
    char buffer[2];
    buffer[1]='\0';
    for (auto c : collected_base) {
        buffer[0]=c;
        *(bIt++) = Rcpp::String(buffer);        
    }
    output[3]=bases;

    return output;
    END_RCPP
}
