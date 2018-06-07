#include "sarlacc.h"

struct rle_walker {
    rle_walker(const char * p, const size_t l) : ptr(p), len(l), 
            last_pos(0), cur_pos(0), nonbases(0), last_base(0), next_base(0),
            true_last_pos(0) {
        while (cur_pos < len) {
            next_base=DNAdecode(ptr[cur_pos]);
            if (next_base!='-') {
                break;
            }
            ++nonbases;
            ++cur_pos;
        }
        return;
    }

    // Advances to the next run length encoding, returning 'false' if we've reached the end.
    void advance() {
        last_pos=cur_pos;
        true_last_pos=last_pos - nonbases; // before 'nonbases' is modified.
        last_base=next_base;
        ++cur_pos;

        while (cur_pos < len) {
            next_base=DNAdecode(ptr[cur_pos]);
            if (next_base!='-' && next_base!=last_base) {
                break;
            }
            ++cur_pos;
            if (next_base=='-') {
                ++nonbases;
            } 
        }
        return;
    }

    bool is_finished () const { return (cur_pos==len); }

    size_t get_start() const { return true_last_pos; }
    size_t get_length()  const { return (cur_pos - nonbases) - true_last_pos; }
    char get_base() const { return last_base; }

    size_t get_start0() const { return last_pos; }
    size_t get_end0() const { return cur_pos; }
private:
    const char * ptr;
    const size_t len;
    size_t last_pos, cur_pos, nonbases, true_last_pos;
    char last_base, next_base;
};

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

        rle_walker SQ(sstr, slen);
        while (!SQ.is_finished()) {
            SQ.advance();
            const size_t homolen=SQ.get_length();
            if (homolen==1) { 
                continue;
            }

            collected_index.push_back(i);
            collected_pos.push_back(SQ.get_start()+1);
            collected_size.push_back(homolen);
            collected_base.push_back(SQ.get_base());
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
        rle_walker SQ(refstr, reflen);
        while (!SQ.is_finished()) {
            SQ.advance();
            const size_t homolen=SQ.get_length();
            if (homolen==1) { 
                continue;
            }

            // Running across everything in the read sequence.
            collected_index.push_back(i);
            collected_pos.push_back(SQ.get_start()+1);
            const size_t left=SQ.get_start0(), right=SQ.get_end0();
            const char curbase=SQ.get_base();
            
            int counter=0;
            for (size_t rpos=left; rpos<right; ++rpos) {
                if (DNAdecode(readstr[rpos])==curbase) {
                    ++counter;
                }
            }
            collected_rlen.push_back(counter);
        }
    }

    return Rcpp::List::create(Rcpp::IntegerVector(collected_index.begin(), collected_index.end()),
        Rcpp::IntegerVector(collected_pos.begin(), collected_pos.end()),
        Rcpp::IntegerVector(collected_rlen.begin(), collected_rlen.end()));
    END_RCPP
}
