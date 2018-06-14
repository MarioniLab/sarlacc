#include "sarlacc.h"
#include "DNA_input.h"

/* Takes local-global alignments (local for the read, global for the adaptor)
 * and restores the characters in the adaptor alignment, while adding gap characters to the read alignment.
 */

SEXP get_aligned_sequence(SEXP aligned_adaptor, SEXP unaligned_adaptor, SEXP start_adaptor, SEXP end_adaptor, SEXP aligned_read) {
    BEGIN_RCPP
    auto adapt_aln=process_DNA_input(aligned_adaptor);
    Rcpp::IntegerVector adapt_starts(start_adaptor), adapt_ends(end_adaptor);

    auto read_aln=process_DNA_input(aligned_read);
    size_t naln=adapt_aln->size();
    if (naln!=read_aln->size()) {
        throw std::runtime_error("number of read and adaptor alignment sequences are not equal");
    }
    
    auto adapt_unaln=process_DNA_input(unaligned_adaptor);
    if (adapt_unaln->size()!=1) {
        throw std::runtime_error("only one unaligned adaptor sequence should be present");
    }
    adapt_unaln->choose(0);
    const char * uptr=adapt_unaln->cstring();
    const size_t ulen=adapt_unaln->length();

    // Creating output objects (guessing the max alignment width, but expanding otherwise).
    Rcpp::StringVector out_adaptors(naln);
    const size_t max_aln=get_max_width(adapt_aln.get());
    std::vector<char> adapt_buffer(max_aln + 1);
    auto set_adaptor =[&] (size_t i, char val) -> void {
        if (i >= adapt_buffer.size()) {
            adapt_buffer.push_back(val);
        } else {
            adapt_buffer[i]=val;
        }
        return;
    };

    Rcpp::StringVector out_reads(naln);
    std::vector<char> read_buffer(max_aln + 1);
    auto set_read = [&] (size_t i, char val) -> void {
        if (i >= read_buffer.size()) {
            read_buffer.push_back(val);
        } else {
            read_buffer[i]=val;
        }
        return;
    };

    // Running through each adaptor string.
    for (size_t a=0; a<naln; ++a) {
        adapt_aln->choose(a);
        read_aln->choose(a);

        const char * aptr=adapt_aln->cstring();
        const size_t alen=adapt_aln->length();

        const char * rptr=read_aln->cstring();
        const size_t rlen=read_aln->length();
        if (rlen!=alen) {
            throw std::runtime_error("read and adaptor alignment strings are not the same length");
        }
       
        // Recreating the full adaptor (global) and read sequence (local, with gaps).
        size_t counter=0;
        const size_t left=adapt_starts[a] - 1;
        while (counter < left) {
            set_adaptor(counter, adapt_unaln->decode(uptr[counter]));
            set_read(counter, '-');
            ++counter;
        }

        for (size_t i=0; i<alen; ++i, ++counter) {
            set_adaptor(counter, adapt_aln->decode(aptr[i]));
            set_read(counter, read_aln->decode(rptr[i]));
        }

        for (size_t i=adapt_ends[a]; i<ulen; ++i, ++counter) {
            set_adaptor(counter, adapt_unaln->decode(uptr[i]));
            set_read(counter, '-');
        }

        set_adaptor(counter, '\0');
        set_read(counter, '\0');
        out_adaptors[a]=Rcpp::String(adapt_buffer.data());
        out_reads[a]=Rcpp::String(read_buffer.data());
    }

    return Rcpp::List::create(out_adaptors, out_reads);
    END_RCPP
}
