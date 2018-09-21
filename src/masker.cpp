#include "masker.h"

masker::masker(double T, Rcpp::NumericVector encoding) : offset(0) {
    Rcpp::StringVector topbases=encoding.names();
    if (topbases.size()==0 || topbases.size()!=encoding.size()) {
        throw std::runtime_error("encoding vector must be non-empty and named");
    }

    char last=0;
    for (size_t i=0; i<topbases.size(); ++i) {
        auto current=Rcpp::as<std::string>(topbases[i]);
        if (current.size()!=1) {
            throw std::runtime_error("names of encoding vector must be one character in length");
        }

        const char curval=current[0];
        if (i > 0) {
            if (curval != last + 1) {
                throw std::runtime_error("names of encoding vector should increase consecutively");
            } else if (encoding[i] > encoding[i-1]) {
                throw std::runtime_error("error probabilities should decrease");
            }
        }

        last=curval;
        if (encoding[i] > T) {
            offset=curval;
        }
    }
    return;
}

void masker::mask(size_t len, const char* seq, const char* qual, char* output) {
    for (size_t counter=0; counter<len; ++counter, ++seq, ++qual) {
        output[counter]=(*qual > offset ? 'N' : *seq);
    }
    return;
}
