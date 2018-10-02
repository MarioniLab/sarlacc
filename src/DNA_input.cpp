#include "DNA_input.h"

DNA_input::DNA_input() : active(NULL), used(0) {}

void DNA_input::clear() {
    used=0;
    return;
}

// String input:

string_input::string_input ( Rcpp::RObject incoming ) : all_values(incoming), holder(1) {};

size_t string_input::size() const {
    return all_values.size();
}

std::pair<const char*, size_t> string_input::get(size_t i) {
    holder[0]=Rcpp::String(all_values[i]);
    const auto& active_string=holder[0];
    used=1;
    return std::make_pair(active_string.get_cstring(), Rf_length(active_string.get_sexp()));
}

std::pair<const char*, size_t> string_input::get_persistent(size_t i) {
    if (used==holder.size()) {
        holder.push_back(Rcpp::String(all_values[i]));
    } else {
        holder[used]=Rcpp::String(all_values[i]);
    }
    const auto& active_string=holder[(used++)];
    return std::make_pair(active_string.get_cstring(), Rf_length(active_string.get_sexp()));
}

size_t string_input::get_len(size_t i) const {
    return Rf_length(Rcpp::String(all_values[i]).get_sexp());
}

// DNAStringSet input:

DNAStringSet_input::DNAStringSet_input ( Rcpp::RObject incoming ) : all_values(hold_XStringSet(SEXP(incoming))), holder(1) {}

size_t DNAStringSet_input::size() const {
    return get_length_from_XStringSet_holder(&all_values);
}

std::pair<const char*, size_t> DNAStringSet_input::get(size_t i) {
    used=1;
    holder[0]=get_internal(0);
    return std::make_pair(holder[0].c_str(), holder[0].size());
}

std::pair<const char*, size_t> DNAStringSet_input::get_persistent(size_t i) {
    if (used==holder.size()) {
        holder.push_back(get_internal(i));
    } else {
        holder[used]=get_internal(i);
    }

    auto& current=holder[used++];
    return std::make_pair(current.c_str(), current.size());
}

std::string DNAStringSet_input::get_internal(size_t i) {
    auto active_string=get_elt_from_XStringSet_holder(&all_values, i);
    const size_t len=active_string.length;
    const char* ptr=active_string.ptr;

    std::string current(len, 0);
    for (size_t i=0; i<len; ++i) {
        current[i]=DNAdecode(ptr[i]);
    }    

    return current;
}

size_t DNAStringSet_input::get_len(size_t i) const {
    auto tmp=get_elt_from_XStringSet_holder(&all_values, i);
    return tmp.length;
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
		const size_t curlen=aln->get_len(i);
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
        const size_t len=seq->get_len(i);
        if (len > ref) {
            ref=len;
        }
    }
    return ref;
}
