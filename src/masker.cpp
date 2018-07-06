#include "masker.h"

masker::masker(double T, size_t mxlen) : buffer(mxlen + 1), threshold(T) {}

Rcpp::String masker::mask(const char* seq, size_t len, double* qual) {
    for (size_t counter=0; counter<len; ++counter, ++seq) {
        buffer[counter]=(qual[counter] > threshold? 'N' : *seq);
    }
    
    buffer[len]='\0';
    return Rcpp::String(buffer.data());
}

unmasker::unmasker(size_t mxlen) : buffer(mxlen + 1) {}

Rcpp::String unmasker::unmask(const char* msk, size_t msklen, const char * origin, size_t orilen) {
	size_t pos_nominal=0;
	for (size_t pos=0; pos<msklen; ++pos) {
		char& outbase=(buffer[pos]=msk[pos]);

		if (outbase!='-') { 
			if (outbase=='N' || outbase=='n') {
				if (pos_nominal >= orilen) {
					throw std::runtime_error("sequence in alignment string is longer than the original");
				}
				outbase=origin[pos_nominal];
			}
			++pos_nominal;
		}
	}

	if (pos_nominal!=orilen) {
		throw std::runtime_error("original sequence and that in the alignment string have different lengths");
		  
	}

	buffer[msklen]='\0';
	return Rcpp::String(buffer.data());
}
