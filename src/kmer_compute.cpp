#include "sarlacc.h"

/**********************************************************************
 * A functor to calculate the k-mer distance between DNAStrings.
 **********************************************************************/

struct kmerize {
private:
    const int klen;
    const XStringSet_holder * ptr;

    int update_kmer (int kmer, const char* seq, int index) {
        if (index >= klen) {
            // Clearing the bits corresponding to the first base in this kmer.
            kmer &= ~(3u << klen*2-2);
        }
        kmer *= 4;
        switch (DNAdecode(seq[index])) {
            case 'A': break; // N's and other ambiguous codes just go here.
            case 'C': kmer += 1; break;
            case 'G': kmer += 2; break;
            case 'T': kmer += 3; break;
        }
        return kmer;
    }

public:
    kmerize(const int k, const XStringSet_holder* p) : klen(k), ptr(p), counts(std::pow(4, klen)) {};

    void process(const int idex) {
        auto sdata=get_elt_from_XStringSet_holder(ptr, idex);
        const char * sptr = sdata.ptr;
        const size_t slen = sdata.length;

        std::fill(counts.begin(), counts.end(), 0);
        int kmer=0, i=0;
        while (i < slen) {
            kmer=update_kmer(kmer, sptr, i);
            ++i;
            if (i < klen && i < slen) { // full k-mers to come.
                continue;
            }
            ++counts[kmer];
        }

        return;
    }

    Rcpp::IntegerVector counts;
};

/**********************************************************************
 * A function to create the k-mer matrix based on their frequencies,
 * and immediately multiplying by another matrix for randomized PCA.
 **********************************************************************/

SEXP get_kmer_matrix (SEXP sequences, SEXP klen) {
    BEGIN_RCPP
    auto seq=hold_XStringSet(sequences); 
    const size_t nseq=get_length_from_XStringSet_holder(&seq);

    const int K=check_integer_scalar(klen, "k-mer length");
    kmerize kmer_counter(K, &seq);
    
    // Setting up the output matrix.
    Rcpp::IntegerMatrix output(nseq, kmer_counter.counts.size());
    auto oIt=output.begin();

    for (size_t i=0; i<nseq; ++i) {
        kmer_counter.process(i);
        auto row=output.row(i);
        std::copy(kmer_counter.counts.begin(), kmer_counter.counts.end(), row.begin());
    }

    return output;
    END_RCPP
}
