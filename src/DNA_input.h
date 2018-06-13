#include "sarlacc.h"
#include <memory>

class DNA_input {
public:    
    DNA_input();
    virtual ~DNA_input();
    const char * cstring() const;

    virtual size_t size() const=0;
    virtual void choose(size_t)=0; 
    virtual size_t length() const=0;

    virtual char decode(char) const=0;
protected:
    const char * active;
};

class string_input : public DNA_input {
public:
    string_input(Rcpp::RObject);
    ~string_input();

    size_t size() const;
    void choose(size_t);
    size_t length() const;
    char decode(char) const;
private:
    Rcpp::StringVector all_values;
    Rcpp::String active_string;
};

class DNAStringSet_input : public DNA_input {
public:
    DNAStringSet_input(Rcpp::RObject);
    ~DNAStringSet_input();

    size_t size() const;
    void choose(size_t);
    size_t length() const;
    char decode(char) const;
private:
    XStringSet_holder all_values;
    Chars_holder active_string;
};

std::unique_ptr<DNA_input> process_DNA_input (Rcpp::RObject);

size_t check_alignment_width(DNA_input * aln);
