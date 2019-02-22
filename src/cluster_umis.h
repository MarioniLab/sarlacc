#ifndef CLUSTER_UMIS_H
#define CLUSTER_UMIS_H

#include "Rcpp.h"

#include <deque>

typedef std::deque<std::deque<int> > value_store;

Rcpp::List cluster_umis(const value_store&);

#endif
