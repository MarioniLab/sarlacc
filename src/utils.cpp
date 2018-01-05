#include "sarlacc.h"

template<typename T, class V> 
T check_scalar(Rcpp::RObject incoming, const char* arg, const char* val) {
    V vec(incoming);
    if (vec.size()!=1) {
        std::stringstream err;
        err << arg << " should be " << val;
        throw std::runtime_error(err.str());
    }
    return vec[0];
}

int check_integer_scalar(Rcpp::RObject incoming, const char* arg) {
    return check_scalar<int, Rcpp::IntegerVector>(incoming, arg, "an integer scalar");
}

double check_numeric_scalar(Rcpp::RObject incoming, const char* arg) {
    return check_scalar<double, Rcpp::NumericVector>(incoming, arg, "a numeric scalar");
}

bool check_logical_scalar(Rcpp::RObject incoming, const char* arg) {
    return check_scalar<bool, Rcpp::LogicalVector>(incoming, arg, "a logical scalar");
}

std::string check_string(Rcpp::RObject incoming, const char* arg) {
    Rcpp::StringVector stuff(incoming);
    check_scalar<Rcpp::String, Rcpp::StringVector>(stuff, arg, "a string");
    return Rcpp::as<std::string>(stuff[0]);
}

int check_alignment_width(XStringSet_holder* aln) {
    const size_t naligns=get_length_from_XStringSet_holder(aln);
    int ref=0;

    for (size_t i=0; i<naligns; ++i) { 
        auto curaln=get_elt_from_XStringSet_holder(aln, i);
        if (i==0) { 
            ref=curaln.length;
        } else if (curaln.length!=ref) {
            throw std::runtime_error("alignment strings should have the same length");
        }
    }
    return ref;
}
