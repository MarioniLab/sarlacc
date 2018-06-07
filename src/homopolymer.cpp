#include "sarlacc.h"

/* This function identifies homopolymers in a given set of sequences,
 * returning the position, length and base of the homopolymer.
 */

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
            continue;
        }

        size_t last_pos=0, cur_pos=0, nonbases=0;
        char last_base=DNAdecode(sstr[0]);

        while ((++cur_pos) < slen) {
            char outbase=DNAdecode(sstr[cur_pos]);

            if (outbase=='-') {
                ++nonbases;
                continue;
            }

            if (outbase==last_base) {
               size_t true_last_pos=last_pos - nonbases; // before 'nonbases' is modified.

               while ((++cur_pos) < slen) {
                    outbase=DNAdecode(sstr[cur_pos]);
                    if (outbase!=last_base) { 
                        break;
                    }
                }

                collected_index.push_back(i);
                collected_pos.push_back(true_last_pos+1);
                collected_size.push_back((cur_pos-nonbases)-true_last_pos);
                collected_base.push_back(last_base);
            }

            // This is always the first element where out_base!=last_base.
            last_pos=cur_pos;
            last_base=outbase;
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

/* This function identifies the homopolymer in a reference alignment string,
 * and returns the number of observed homopolymers in the read alignment string.
 * Some effort is required to handle the deletion character '-'.
 */

SEXP match_homopolymers (SEXP ref_align, SEXP read_align) {
    BEGIN_RCPP

    auto rf_aln=hold_XStringSet(ref_align); 
    const size_t nalign=get_length_from_XStringSet_holder(&rf_aln);
    auto rd_aln=hold_XStringSet(read_align); 
    if (nalign!=get_length_from_XStringSet_holder(&rd_aln)) {
        throw std::runtime_error("lengths of alignment vectors should match up");
    }
  
    // Setting up output containers.
    std::deque<int> collected_index, collected_pos, collected_rlen;

    // Running through the pairwise aligners and reporting the observed length of all homopolymers.
    for (size_t i=0; i<nalign; ++i) {
        auto curref=get_elt_from_XStringSet_holder(&rf_aln, i);
        const char* refstr=curref.ptr;
        const size_t reflen=curref.length;

        auto curread=get_elt_from_XStringSet_holder(&rd_aln, i);
        const char* readstr=curread.ptr;
        const size_t readlen=curread.length;

        if (readlen!=reflen) {
            throw std::runtime_error("read and reference alignment strings should have equal length");
        }
        if (reflen==0) {
            continue;
        }

        // Iterating across the reference and counting the number of homopolymers in the read.
        size_t last_pos=0, cur_pos=0, nonbases=0;
        char last_base=DNAdecode(refstr[0]); // doesn't matter that it's a '-', as it will get overwritten in the loop.

        while ((++cur_pos) < reflen) {
            char outbase=DNAdecode(refstr[cur_pos]);

            if (outbase=='-') {
                ++nonbases;
                continue;
            }

            if (outbase==last_base) {
                size_t true_last_pos=last_pos - nonbases; // before 'nonbases' is modified.

                while ((++cur_pos) < reflen) {
                    outbase=DNAdecode(refstr[cur_pos]);
                    if (outbase=='-') { 
                        ++nonbases;
                        continue;
                    } 
                    if (outbase!=last_base) { 
                        break;
                    }
                }

                // Running across everything in the read sequence.
                collected_index.push_back(i);
                collected_pos.push_back(true_last_pos+1);
                int counter=0;
                for (size_t rpos=last_pos; rpos<cur_pos; ++rpos) {
                    if (DNAdecode(readstr[rpos])==last_base) {
                        ++counter;
                    }
                }
                collected_rlen.push_back(counter);
            }

            // This is always the first element where out_base!=last_base.
            last_pos=cur_pos;
            last_base=outbase;
        }
    }

    return Rcpp::List::create(Rcpp::IntegerVector(collected_index.begin(), collected_index.end()),
        Rcpp::IntegerVector(collected_pos.begin(), collected_pos.end()),
        Rcpp::IntegerVector(collected_rlen.begin(), collected_rlen.end()));
    END_RCPP
}
