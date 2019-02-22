#ifndef UMI_CLUSTERER_H
#define UMI_CLUSTERER_H

#include "Rcpp.h"
#include "value_store.h"
#include <deque>

class umi_clusterer{
public:
    umi_clusterer() = default;
    value_store<int, std::deque<int> > storage;
    Rcpp::List cluster();
private:
    std::deque<size_t> ordering, remaining;
    std::deque<int> per_cluster_values;
    std::deque<Rcpp::IntegerVector> output;
};

#endif
