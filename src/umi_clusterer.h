#ifndef UMI_CLUSTERER_H
#define UMI_CLUSTERER_H
#include "sarlacc.h"

class umi_clusterer{
public:
    umi_clusterer();

    template<class IT>
    void add(IT start, IT end) {
        if (n_stored > storage_start.size()) {
            storage_start.resize(n_stored);
            storage_len.resize(n_stored);
        }

        const auto len=end - start;
        storage_start[n_stored]=nvals_stored;
        storage_len[n_stored]=len;

        const auto required=nvals_stored + len;
        if (required > storage_values.size()) {
            storage_values.resize(required);
        }

        std::copy(start, end, storage_values.begin() + nvals_stored);
        nvals_stored=required;
        ++n_stored;
        return;
    }

    void clear();

    Rcpp::List cluster();
private:
    size_t n_stored, nvals_stored;
    std::deque<int> storage_values, storage_start, storage_len;

    std::deque<size_t> ordering, remaining;

    std::deque<int> per_cluster_values;
    std::deque<Rcpp::IntegerVector> output;
};

#endif
