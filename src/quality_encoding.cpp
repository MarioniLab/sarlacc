#include "quality_encoding.h"

quality_encoding::quality_encoding(Rcpp::NumericVector encoding) : offset(0), errors(encoding) {
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
        } else {
            offset=curval;
        }

        last=curval;
    }
    return;
}

char quality_encoding::lowest() const { return offset; }

Rcpp::NumericVector quality_encoding::get_errors() const { return errors; }

double quality_encoding::to_error(char Q) const {
    if (Q < offset) {
        throw std::runtime_error("quality cannot be lower than smallest encoded value");
    }
    size_t i=Q - offset;
    if (i > errors.size()) {
        i=errors.size() - 1; // must be non-empty, so this won't underflow.
    }
    return errors[i];
}
