#include "sarlacc.h"

/* This function identifies errors at each base in a reference alignment string. */

SEXP find_errors (SEXP ref_align, SEXP read_align) {
    BEGIN_RCPP

    auto rf_aln=hold_XStringSet(ref_align); 
    const size_t nalign=get_length_from_XStringSet_holder(&rf_aln);
    auto rd_aln=hold_XStringSet(read_align); 
    if (nalign!=get_length_from_XStringSet_holder(&rd_aln)) {
        throw std::runtime_error("lengths of alignment vectors should match up");
    }

    // Setting up output containers.
    size_t standard_len=0;
    if (nalign) {
        auto standard=get_elt_from_XStringSet_holder(&rf_aln, 0);
        const char* sstr=standard.ptr;
        const size_t slen=standard.length;
        for (size_t s=0; s<slen; ++s) {
            if (DNAdecode(sstr[s])!='-') { ++standard_len; }
        }
    }
    Rcpp::IntegerVector all_deletions(standard_len), all_to_A(standard_len), all_to_C(standard_len), 
        all_to_G(standard_len), all_to_T(standard_len);
    std::deque<int> insertion_pos, insertion_len;

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

        // Iterating across the reference.
        size_t cur_pos=0, nonbases=0;
        while (cur_pos < reflen) {
            const char ref_base=DNAdecode(refstr[cur_pos]), read_base=DNAdecode(readstr[cur_pos]);

            if (ref_base!=read_base) {
                if (ref_base!='-') {
                    // Adding to the statistics directly, if we're not currently in an insertion.
                    const size_t true_pos=cur_pos - nonbases;
                    if (true_pos >= standard_len) {
                        throw std::runtime_error("reference sequence should be the same for all alignments");
                    }

                    switch (read_base) {
                        case '-':
                            ++(all_deletions[true_pos]);
                            break;
                        case 'A':
                            ++(all_to_A[true_pos]);
                            break;
                        case 'C':
                            ++(all_to_C[true_pos]);
                            break;
                        case 'G':
                            ++(all_to_G[true_pos]);
                            break;
                        case 'T':
                            ++(all_to_T[true_pos]);
                            break;
                        default:
                            std::stringstream err;
                            err << "unknown character '" << read_base << "' in alignment string";
                            throw std::runtime_error(err.str().c_str());
                    }
                    ++cur_pos;
                } else {
                    // Determining the size of the insertion and adding it at the position of the next base.
                    // This may be at the end of the sequence. 
                    const size_t previous=cur_pos;
                    ++cur_pos;
                    ++nonbases;

                    while (cur_pos < reflen) {
                        const char cur_base=DNAdecode(refstr[cur_pos]);
                        if (cur_base!='-') { 
                            break;
                        }
                        ++cur_pos;
                        ++nonbases;
                    }

                    const size_t true_pos=cur_pos - nonbases;
                    insertion_pos.push_back(true_pos);
                    insertion_len.push_back(cur_pos - previous);
                }
            } else {
                ++cur_pos;
                if (ref_base=='-') {
                    ++nonbases;
                }
            }
        }
    }

    return Rcpp::List::create(all_to_A, all_to_C, all_to_G, all_to_T, all_deletions, 
        Rcpp::IntegerVector(insertion_pos.begin(), insertion_pos.end()),
        Rcpp::IntegerVector(insertion_len.begin(), insertion_len.end()));
    END_RCPP
}
