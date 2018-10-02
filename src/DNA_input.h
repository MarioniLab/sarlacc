#ifndef DNA_INPUT_H
#define DNA_INPUT_H

#include "sarlacc.h"
#include <memory>
#include <string>

class DNA_input {
public:    
    DNA_input();
    virtual ~DNA_input() = default;
    
    virtual size_t size() const=0;
    virtual std::pair<const char*, size_t> get(size_t)=0; 
    virtual std::pair<const char*, size_t> get_persistent(size_t)=0; 
    virtual size_t get_len(size_t) const=0;
    
    virtual void clear();
protected:
    const char * active;
    size_t used;
};

class string_input : public DNA_input {
public:
    string_input(Rcpp::RObject);
    ~string_input() = default;

    size_t size() const;
    std::pair<const char*, size_t> get(size_t); 
    std::pair<const char*, size_t> get_persistent(size_t); 
    size_t get_len(size_t) const;
private:
    Rcpp::StringVector all_values;
    std::deque<Rcpp::String> holder;
};

class DNAStringSet_input : public DNA_input {
public:
    DNAStringSet_input(Rcpp::RObject);
    ~DNAStringSet_input() = default;

    size_t size() const;
    std::pair<const char*, size_t> get(size_t); 
    std::pair<const char*, size_t> get_persistent(size_t); 
    size_t get_len(size_t) const;
private:
    XStringSet_holder all_values;
    std::deque<std::string> holder;
    std::string get_internal(size_t);
};

std::unique_ptr<DNA_input> process_DNA_input (Rcpp::RObject);

size_t check_alignment_width(DNA_input * aln);

size_t get_max_width(DNA_input * aln);

#endif
