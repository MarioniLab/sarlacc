#ifndef DNA_INPUT_H
#define DNA_INPUT_H

#include "sarlacc.h"
#include <memory>

class DNA_input {
public:    
    DNA_input();
    virtual ~DNA_input();
    
    virtual size_t size() const=0;
    virtual std::pair<const char*, size_t> get(size_t)=0; 
    virtual size_t get_len(size_t) const=0;
    virtual void clear()=0;
protected:
    const char * active;
};

class string_input : public DNA_input {
public:
    string_input(Rcpp::RObject);
    ~string_input();

    size_t size() const;
    std::pair<const char*, size_t> get(size_t); 
    size_t get_len(size_t) const;
    void clear();
private:
    Rcpp::StringVector all_values;
    std::deque<Rcpp::String> holder;
};

class DNAStringSet_input : public DNA_input {
public:
    DNAStringSet_input(Rcpp::RObject);
    ~DNAStringSet_input();

    size_t size() const;
    std::pair<const char*, size_t> get(size_t); 
    size_t get_len(size_t) const;
    void clear();
private:
    XStringSet_holder all_values;
    std::deque<Chars_holder> holder;
    std::vector<char> buffer;
    size_t used;
};

std::unique_ptr<DNA_input> process_DNA_input (Rcpp::RObject);

size_t check_alignment_width(DNA_input * aln);

size_t get_max_width(DNA_input * aln);

#endif
