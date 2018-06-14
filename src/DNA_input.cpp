#include "DNA_input.h"

DNA_input::DNA_input() : active(NULL) {}

DNA_input::~DNA_input() {}

const char * DNA_input::cstring() const { return active; }

string_input::string_input ( Rcpp::RObject incoming ) : all_values(incoming) {};

string_input::~string_input() {}

size_t string_input::size() const {
    return all_values.size();
}

/* You need to "choose" the entry before you can access the C-string pointer.
 * This makes a copy of the intermediate object, which means that the pointer 
 * must be valid as long as the intermediate object is still alive.
 */ 
void string_input::choose(size_t i) {
    active_string=all_values[i];
    active=active_string.get_cstring();
    return;
}

size_t string_input::length() const {
    return Rf_length(active_string.get_sexp());
}

char string_input::decode(char incoming) const {
    return incoming;
}

DNAStringSet_input::DNAStringSet_input ( Rcpp::RObject incoming ) : all_values(hold_XStringSet(SEXP(incoming))) {}

DNAStringSet_input::~DNAStringSet_input() {}

size_t DNAStringSet_input::size() const {
    return get_length_from_XStringSet_holder(&all_values);
}

void DNAStringSet_input::choose(size_t i) {
    active_string=get_elt_from_XStringSet_holder(&all_values, i);
    active=active_string.ptr;
    return;
}

size_t DNAStringSet_input::length() const {
    return active_string.length;
}

char DNAStringSet_input::decode(char incoming) const {
    return DNAdecode(incoming);
}

std::unique_ptr<DNA_input> process_DNA_input (Rcpp::RObject incoming) {
    if (incoming.isS4()) {
        return std::unique_ptr<DNA_input>(new DNAStringSet_input(incoming));
    } else {
        return std::unique_ptr<DNA_input>(new string_input(incoming));
    }
}

size_t check_alignment_width(DNA_input * aln) {
    const size_t naligns=aln->size();
    size_t ref=0;

    for (size_t i=0; i<naligns; ++i) { 
		aln->choose(i);
		const size_t curlen=aln->length();

        if (i==0) { 
            ref=curlen;
        } else if (curlen!=ref) {
            throw std::runtime_error("alignment strings should have the same length");
        }
    }

    return ref;
}

size_t get_max_width(DNA_input * seq) {
    const size_t nseq=seq->size();
    size_t ref=0;
    for (size_t i=0; i<nseq; ++i) { 
        seq->choose(i);
        const size_t len=seq->length();
        if (len > ref) {
            ref=len;
        }
    }
    return ref;
}
